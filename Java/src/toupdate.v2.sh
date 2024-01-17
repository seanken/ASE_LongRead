cp *.java ../Gradle/app/src/main/java/allelemasseq
cd ../G*
./gradlew build
cp app/build/libs/AlleleMASSeq.jar ../../scripts/AlleleMASSeq.jar
cd -
