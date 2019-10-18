package org.chusj;

import redis.clients.jedis.Jedis;
import redis.clients.jedis.JedisPool;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class RedisGeneSetHelper {

//    private static Jedis jedisClient = new Jedis("localhost", 6379);
    private static final String redisHost = "localhost";
    private static final Integer redisPort = 6379;
    private static final int NUMBER_OF_RETRY = 10;
    private static Map<String,Set<String>> redisCache = new HashMap<>();
    private static JedisPool pool = new JedisPool(redisHost, redisPort);

    public static void main(String[] args) {

        Set<String> geneSets = getMembersForEnsId("ENSG00000142611");

        System.out.println("Member for "+ "ENSG00000142611 are:" );
        geneSets.forEach( (String member) -> System.out.println("\t"+member ));

        Set<String> geneSetsFromCache = getMembersForEnsId("ENSG00000142611");

        System.out.println("Member from cache? for "+ "ENSG00000142611 are:" );
        geneSetsFromCache.forEach( (String member) -> System.out.println("\t"+member ));

    }

    static Set<String> getMembersForEnsId(String ensId) {

        Set<String> sMembers = new HashSet<>();
        if (redisCache.containsKey(ensId)) {
          return redisCache.get(ensId);
        }
        for (int i=0; i < NUMBER_OF_RETRY; i++)
            try (Jedis jedisClient = pool.getResource()) {
                sMembers = jedisClient.smembers("id:" + ensId);
                redisCache.put(ensId,sMembers);
                break;
            } catch (Exception e) {
                if (i< (NUMBER_OF_RETRY - 1) ) {
                    //System.err.println("Redis cle #" + i + " failed... Retrying...");
                try {
                    Thread.sleep(100);
                } catch (InterruptedException e1) {
                    e1.printStackTrace();
                }
                } else {
                    System.err.println("XXXXX Redis call failed "+NUMBER_OF_RETRY + " times");
                    e.printStackTrace();
                }
            }
        return sMembers;

    }

    static Set<String> getIdsForGeneSymbol(String geneSymbol) {

        Set<String> sMembers = new HashSet<>();
        for (int i=0; i< 3; i++) {
            try (Jedis jedisClient = pool.getResource()){
                sMembers = jedisClient.smembers("gene:"+geneSymbol);
                break;
            } catch (Exception e) {
                System.err.println("Redis call #"+i+" failed... Retrying...");
            }
        }

        return sMembers;

    }


}
