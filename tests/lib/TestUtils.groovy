import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.StandardCopyOption

class TestUtils {
    // to be able to refer to the cache dir in setup blocks
    static String getCacheDir(Object baseDir) {
        return "${baseDir}/.nf-test/shared_cache"
    }

    static void copyFolder(String source_string, String destination_string) {
        def source = Paths.get(source_string)
        def destination = Paths.get(destination_string)
        Files.walk(source).each { src ->
            Path dest = destination.resolve(source.relativize(src))
            if (Files.isDirectory(src)) {
                Files.createDirectories(dest)
            } else {
                Files.copy(src, dest, StandardCopyOption.REPLACE_EXISTING)
            }
        }
    }

    static void copyFiles(Map<String, String> files_md5, String destination_string) {
        print "\n  INFO: Copying files to ${destination_string}..."
        def destination = Paths.get(destination_string)
        Files.createDirectories(destination)
        files_md5.each { filepath, expected_md5 ->
            def src = Paths.get(filepath)
            def dest = destination.resolve(src.getFileName())
            Files.copy(src, dest, StandardCopyOption.REPLACE_EXISTING)
            def actual_md5 = "md5sum ${dest}".execute().text.trim().split()[0]
            if (actual_md5 != expected_md5) {
                throw new RuntimeException("MD5 mismatch for ${filepath}: expected ${expected_md5}, got ${actual_md5}")
            }
        }
        print "\n  INFO: Files copied to ${destination_string}"
    }

    static void checkMD5s(Map<String, String> md5s) {
        if ( !md5s ) {
            throw new Exception( "md5s not provided. Cannot verify integrity of cached archive." )
        }
        md5s.each { file, expected_md5 ->
            if ( new File(file).exists() ) {
                def actual_md5 = "md5sum ${file}".execute().text.trim().split()[0]
                if ( actual_md5 != expected_md5 ) {
                    throw new Exception( "md5 mismatch for ${file}. Expected ${expected_md5}, got ${actual_md5}. Please clean the cache and retry." )  
                }
            } else {
                throw new Exception( "Expected file ${file} not found in cache. Please clean the cache and retry." )
            }
        }
    }

    static void stageTestData(Object baseDir) {
        def source = "${baseDir}/assets/test-data/data"
        def destination = "${getCacheDir(baseDir)}/data"
        print "\n  INFO: Downloads completed. Copying test data from ${source} to ${destination}..."
        copyFolder( source, destination )
        print "\n  INFO: Test data copied to ${destination}"
    }

    static void fetchArchive (List<String> urls, String destination, String touch_file = null, Map<String, String> md5s) {
        if ( new File( destination ).exists() && new File(destination).isDirectory() ) {
            print "\n  INFO: Using cached archive at ${destination}"
        } else {
            print "\n  INFO: Cache missing. Downloading archive to ${destination}..."
            urls.each { url -> curlAndExtract( url, destination ) }
            if ( touch_file ) {
                def marker = new File( touch_file )
                marker << ""
            }
        }
        checkMD5s( md5s )
        print "\n  INFO: Cached archive at ${destination} passed integrity check."
    }

    static void makeDatabaseSheet(Object launchDir, String content) {
        def database_csv = new File( "${launchDir}/database.csv" )
        database_csv.text = "tool,db_name,db_params,db_path\n" + content
    }

    static void fetchHostremovalT2TGenome(Object baseDir) {
        def hostremoval_t2t_genome = "${getCacheDir(baseDir)}/hostremoval/t2t_chm13v2/"
        def t2t_fa_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz"
        def t2t_fa_destination = new File( "${hostremoval_t2t_genome}/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz" )
        def bowtie2_index_url = "https://genome-idx.s3.amazonaws.com/bt/chm13v2.0.zip"
        if ( new File( hostremoval_t2t_genome ).exists() && new File( hostremoval_t2t_genome ).isDirectory() ) {
            print "\n  INFO: Using cached reference genome at ${hostremoval_t2t_genome}"
        } else {
            print "\n  INFO: Cache missing. Downloading reference genome to ${hostremoval_t2t_genome}..."
            new File( hostremoval_t2t_genome ).mkdirs()
            t2t_fa_destination.withOutputStream { out -> out << t2t_fa_url.toURL().openStream() }
            curlAndExtract( bowtie2_index_url, "${hostremoval_t2t_genome}/bowtie2_index" )
        }
        def md5s = [
            "${hostremoval_t2t_genome}/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz": "9e6bf6b586bc8954208d1cc1d5f2fc99",
            "${hostremoval_t2t_genome}/bowtie2_index/chm13v2.0/chm13v2.0.1.bt2": "20901bcfe2065fc399d1d53ff87d39b3",
            "${hostremoval_t2t_genome}/bowtie2_index/chm13v2.0/chm13v2.0.2.bt2": "5f09ff75d61ecae9a14fbf26916f2bcf",
            "${hostremoval_t2t_genome}/bowtie2_index/chm13v2.0/chm13v2.0.3.bt2": "202208bb6de8fb51010d50a7e56aa717",
            "${hostremoval_t2t_genome}/bowtie2_index/chm13v2.0/chm13v2.0.4.bt2": "798cadba311a10d612b877b2c5382536",
            "${hostremoval_t2t_genome}/bowtie2_index/chm13v2.0/chm13v2.0.rev.1.bt2": "b3edaa2215e7d87d7306105789b126c2",
            "${hostremoval_t2t_genome}/bowtie2_index/chm13v2.0/chm13v2.0.rev.2.bt2": "95f530d31dc402e9af9898874ff99c81",
        ]
        checkMD5s( md5s )
        print "\n  INFO: Reference genome at ${hostremoval_t2t_genome} passed integrity check."
    }

    static void generateBwaIndex(String fasta, Map<String, String> md5s = null) {
        def indexFiles = [".amb", ".ann", ".bwt", ".pac", ".sa"]
        if ( indexFiles.every { new File("${fasta}${it}").exists() } ) {
            print "\n  INFO: Using cached BWA index for ${fasta}."
        } else {
            print "\n  INFO: Generating BWA index for ${fasta}..."
            print "\n  WARN: BWA 0.7.19-r1273 must be installed and in your PATH for this to work!"
            def proc = "bwa index ${fasta}".execute()
            def stderr = new StringBuilder()
            proc.consumeProcessErrorStream(stderr)
            def exitCode = proc.waitFor()
            if ( exitCode != 0 ) {
                throw new Exception( "Failed to generate BWA index for ${fasta}." )
            }
        }
        checkMD5s( md5s )
        print "\n  INFO: BWA index for ${fasta} passed integrity check."
    }

    static List<String> getUnstableName() {
        return [ 'pipeline_info/*.{html,json,txt}' ]
    }
}