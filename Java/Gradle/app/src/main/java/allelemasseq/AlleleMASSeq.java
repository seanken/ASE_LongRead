package AlleleMASSeq;
import java.io.FileWriter;
import java.lang.*;
import java.io.*;
import htsjdk.samtools.*;
import java.util.*;
import java.util.Map.Entry;
import java.util.Map;


public class AlleleMASSeq
{

    public static void main(String args[]) throws Exception
    {
        Long_ASE_Count ASECounter=new Long_ASE_Count(args);
        ASECounter.IterateOverAll();
    }

}

