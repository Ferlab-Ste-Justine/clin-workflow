package org.chusj;

import redis.clients.jedis.Jedis;

import java.util.HashSet;
import java.util.Set;

public class RedisGeneSetHelper {

    private static Jedis jedisClient = new Jedis();
    private static final int NUMBER_OF_RETRY = 3;

    public static void main(String[] args) {

        Set<String> geneSets = getMembersForEnsId("ENSG00000166813");

        System.out.println("Member for "+ "ENSG00000166813 are:" );

        geneSets.forEach( (String member) -> System.out.println("\t"+member ));


    }

    static Set<String> getMembersForEnsId(String ensId) {


        Set<String> sMembers = new HashSet<>();
        for (int i=0; i< NUMBER_OF_RETRY; i++)
            try {
                sMembers = jedisClient.smembers("id:" + ensId);
                break;
            } catch (Exception e) {
                if (i<NUMBER_OF_RETRY) {
                    System.err.println("Redis cle #" + i + " failed... Retrying...");
                try {
                    Thread.sleep(400);
                } catch (InterruptedException e1) {
                    e1.printStackTrace();
                }
                } else {
                    System.err.println("Redis call failed "+NUMBER_OF_RETRY + " times");
                }
            }

        return sMembers;

    }

    static Set<String> getIdsForGeneSymbol(String geneSymbol) {

        Set<String> sMembers = new HashSet<>();
        for (int i=0; i< 3; i++) {
            try {
                sMembers = jedisClient.smembers("gene:"+geneSymbol);
                break;
            } catch (Exception e) {
                System.err.println("Redis call #"+i+" failed... Retrying...");
            }
        }

        return sMembers;

    }


}
