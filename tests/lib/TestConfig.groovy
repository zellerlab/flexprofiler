// to be able to refer to the cache dir in setup blocks

class TestConfig {
    static String getCacheDir(Object baseDir) {
        return "${baseDir}/.nf-test/shared_cache"
    }
}