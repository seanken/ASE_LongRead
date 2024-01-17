Jppackage AlleleMASSeq;
import java.io.FileWriter;
import java.lang.*;
import java.io.*;
import htsjdk.samtools.*;
import java.util.*;
import java.util.Map.Entry;
import java.util.Map;


public class Long_ASE_Count
{

    protected FileWriter savFil;
    protected File bamFile;
    protected HashMap<String, Boolean> map;
    protected HashMap<String, String> UMIMap;
    protected HashMap<String,Integer> counts;
    protected int Good;
    protected int Bad;
    protected int numRem;
    protected String gene_tag;
    protected String umi_tag;

    public Long_ASE_Count(String args[]) throws Exception
    {
        //ArrayList<String[]> ret=new ArrayList<String[]>();
        this.bamFile = new File(args[0]); 
        File vcfFile=new File(args[1]);
        //int ind=Integer.parseInt(args[2]);
        int ind=0; //assums vcf is prepared 
        String saveFile=args[2];
        this.gene_tag=args[3];
        this.umi_tag=args[4];


        this.savFil = new FileWriter(saveFile); 

        print("Create VCF HashMap");

        this.map= HashMapVCF(vcfFile,ind);
        System.out.println(this.map.size());
    }


    public void IterateOverAll() throws Exception
    {

        print("Read in long read data!");
        SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(bamFile);
        SAMRecordIterator r = sr.iterator();
        int iter=0;
        this.Good=0;
        this.Bad=0;
        this.numRem=0;

        this.UMIMap=new HashMap<String, String>(); 

        while(r.hasNext()) {
            iter=iter+1;

            if(iter%10000==0)
            {
                System.out.println(this.Good);
                System.out.println(this.Bad);
                System.out.println(this.numRem);
                System.out.println((float)this.Good/(float)(this.Good+this.Bad));
                System.out.println((float)this.Good/(float)(iter));
                System.out.println(iter);
            }

            SAMRecord read=r.next();

            this.ProcessRead(read);

        }

        System.out.println("Get Allele Counts");
        this.UMIMAPtoCounts();
        
        System.out.println("Save!");
        this.saveCounts();

        
        r.close();
        sr.close();
        this.savFil.close();

    }




    //Saves the UMI counts 
    public void saveCounts() throws Exception
    {
        Iterator<String> it_cnt=this.counts.keySet().iterator();

        while (it_cnt.hasNext()) {
            String key=it_cnt.next();
            Integer val=this.counts.get(key);

            String res=key+" "+Integer.toString(val);
            this.savFil.write(res+ System.lineSeparator());

        }

    }

    //Takes the internal UMIMAP object and gets UMI count information from it
    public void UMIMAPtoCounts()
    {
        this.counts=new HashMap<String,Integer>();
        Iterator<String> it=this.UMIMap.keySet().iterator();


        while (it.hasNext()) {
            String key=it.next();
            String val=this.UMIMap.get(key);
            String[] split=key.split(" ");
            String cbc=split[1];
            String gene=split[2];
            //String allele=split[3];
            String res=cbc+" "+gene+" "+val;
            //savFil.write(res+ System.lineSeparator());
            if(!this.counts.containsKey(res))
            {
                this.counts.put(res,0);
            }

            this.counts.put(res,this.counts.get(res)+1);

        }
    }



    public static HashMap<String,Boolean> HashMapVCF(File vcfFile,int ind) throws Exception
    {

     
        Scanner sc = new Scanner(vcfFile); 

        HashMap<String,Boolean> vcfMap=new HashMap<String,Boolean>();

        int iter=0;

        while (sc.hasNextLine()) 
        {
            iter=iter+1;

            if(iter %100000==0)
            {
                System.out.println(iter);
            }

            String line=sc.nextLine();
            if('#'==line.charAt(0))
            {
                continue;
            }

            String[] split_line=line.split("\\s+");


            String chrom=split_line[0];
            String pos=split_line[1];
            String ref=split_line[3];
            String alt=split_line[4];
            int genoPos=9+ind;
            String geno=split_line[genoPos];
            String[] split_geno=geno.split(":");
            geno=split_geno[0];

            if(ref.length()>1 || alt.length()>1)
            {
                continue;
            }
                

            chrom=chrom.replace("chr","");

            String key1=chrom+"_"+pos+"_"+ref;
            String key2=chrom+"_"+pos+"_"+alt;

            Boolean val1=true;
            Boolean val2=false;

            if(geno.equals("0|1") || geno.equals("0/1"))
            {
                vcfMap.put(key1,val1);
                vcfMap.put(key2,val2);
            }


            if(iter %100000==0)
            {
                System.out.println(key1);
                System.out.println("");
            }




            if(geno.equals("1|0") || geno.equals("1/0"))
            {
                vcfMap.put(key1,val2);
                vcfMap.put(key2,val1);
            }


        }

        return(vcfMap);

    }



    public void ProcessRead(SAMRecord read)
    {
        String seq=read.getReadString();

        String cbc=read.getStringAttribute(this.gene_tag);
        String umi=read.getStringAttribute(this.umi_tag);
        String gene=read.getStringAttribute("XT");
        Character isprim=read.getCharacterAttribute("tp");
        String suppmaps=read.getStringAttribute("SA");

        if(!(suppmaps==null))
        {
            return;
        }

        //Add in multithis.map filtering

        if(gene==null || cbc==null || umi==null || isprim==null)
        {
            return;
        }

        if(isprim!='P')
        {
            return;
        }
        this.numRem=this.numRem+1;

        //getReferencePositionAtReadPosition
        int start=read.getStart();
        String chrom=read.getContig();
        chrom=chrom.replace("chr","");
        int readLen=read.getReadLength();
        int All1=0;
        int All2=0;


        for(int i=0;i<readLen;i=i+1)
        {

            int pos=read.getReferencePositionAtReadPosition(i+1);
            //pos=pos;

            if(pos<1){continue;}


            String base=Character.toString(seq.charAt(i));


            String key=chrom+"_"+Integer.toString(pos)+"_"+base;

            //if(iter % 10000==0)
            //{
            //System.out.println(key);
            //}


            if(this.map.containsKey(key))
            {
                Boolean val=this.map.get(key);
                if(val)
                {
                    All1=All1+1;
                }
                else
                {
                    All2=All2+1;
                }

            }


        }

        int Tot=All1+All2;
        if(Tot<1){return;}
        float rat=(float)Math.max(All1,All2)/(float)Tot;

        if(rat>.95)
        {
            this.Good=this.Good+1;
        }
        else
        {
            this.Bad=this.Bad+1;
            return;
        }



        String res=umi+" "+cbc+" "+gene;
        String val="None";
        if(All1>All2)
        {
            val="All1";
        }
        if(All2>All1)
        {
            val="All2";
        }

        if(this.UMIMap.containsKey(res))
        {
            String cur=this.UMIMap.get(res);
            if(!cur.equals(val))
            {
            this.UMIMap.put(res,"Ambig");
            }
        }
        else
        {
            this.UMIMap.put(res,val);
        }
    }



    public static void print(String toPrint)
    {

    System.out.println(toPrint);

    }

}

