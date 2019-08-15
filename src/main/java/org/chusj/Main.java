package org.chusj;



public class Main {
    public static void main(String[] args) throws Exception {

        if (args.length == 0 ) {
            MemcachedJavaClient.main(args);
        } else if ( args.length == 2) {
            VepHelper.main(args);
        } else {
            Validator.main(args);
        }
    }
}
