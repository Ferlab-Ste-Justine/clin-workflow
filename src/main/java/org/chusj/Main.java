package org.chusj;



public class Main {
    public static void main(String[] args) throws Exception {

        if (args.length ==0 ) {
            MemcachedJavaClient.main(args);
            return;
        } else {
            Validator.main(args);
        }
    }
}
