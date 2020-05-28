package org.chusj;



public class Main {
    public static void main(String[] args) throws Exception {

        if (args.length == 0 ) {
            RedisGeneSetHelper.main(args);
        } else if ( args.length == 2) {
            RedisGeneSetIndexer.main(args);
        } else if ( args.length == 3) {
            VepHelper.main(args);
        } else if (args.length == 4) {
            Validator.main(args);
        }
    }
}
