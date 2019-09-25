package org.chusj;

import redis.clients.jedis.Jedis;

import java.util.HashSet;
import java.util.Set;

public class RedisGeneSetHelper {

    private static Jedis jedisClient = new Jedis();

    public static void main(String[] args) {

        Set<String> geneSets = getMembersForEnsId("ENSG00000166813");

        System.out.println("Member for "+ "ENSG00000166813 are:" );

        geneSets.forEach( (String member) -> System.out.println("\t"+member ));


    }

    static Set<String> getMembersForEnsId(String ensId) {

        Set<String> sMembers = new HashSet<>();
        for (int i=0; i< 3; i++) {
            try {
                sMembers = jedisClient.smembers("id:"+ensId);
                break;
            } catch (Exception e) {
                System.err.println("Redis call #"+i+" failed... Retrying...");
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
