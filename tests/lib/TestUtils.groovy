import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.StandardCopyOption

class TestUtils {
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

    // to be able to refer to the cache dir in setup blocks
    static String getCacheDir(Object baseDir) {
        return "${baseDir}/.nf-test/shared_cache"
    }
}