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

    static void stageTestData(Object baseDir) {
        def source = "${baseDir}/assets/test-data/data"
        def destination = "${getCacheDir(baseDir)}/data"
        print "\n  INFO: Downloads completed. Copying test data from ${source} to ${destination}..."
        copyFolder( source, destination )
        print "\n  INFO: Test data copied to ${destination}"
    }

    static void fetchArchive (String url, String destination, String touch_file = null) {
            if ( new File( destination ).exists() && new File(destination).isDirectory() ) {
                print "\n  INFO: Using cached archive at ${destination}"
            } else {
                print "\n  INFO: Cache missing. Downloading archive to ${destination}..."
                curlAndExtract( url, destination )
                if ( touch_file ) {
                    def marker = new File( touch_file )
                    marker << ""
                }
            }
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
    }

    static String getUnstableName() {
        return "pipeline_info/execution_*.{html,txt}"
    }
}