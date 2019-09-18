package org.chusj;

import redis.clients.jedis.Jedis;

import java.util.Set;

public class RedisGeneSetHelper {

    private static Jedis jedisClient = new Jedis();

    public static void main(String[] args) {

        Set<String> geneSets = getMembersForEnsId("ENSG00000166813");

        System.out.println("Member for "+ "ENSG00000166813 are:" );

        geneSets.forEach( (String member) -> System.out.println("\t"+member ));


    }

    static Set<String> getMembersForEnsId(String ensId) {

        return jedisClient.smembers("id:"+ensId);

    }
}
