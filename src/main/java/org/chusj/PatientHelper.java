package org.chusj;

import org.apache.http.HttpHost;
import org.elasticsearch.action.get.GetRequest;
import org.elasticsearch.action.get.GetResponse;
import org.elasticsearch.action.search.SearchRequest;
import org.elasticsearch.action.search.SearchResponse;
import org.elasticsearch.client.RequestOptions;
import org.elasticsearch.client.RestClient;
import org.elasticsearch.client.RestHighLevelClient;
import org.elasticsearch.index.query.QueryStringQueryBuilder;
import org.elasticsearch.search.builder.SearchSourceBuilder;
import org.json.JSONArray;
import org.json.JSONObject;

import java.util.*;

import static org.elasticsearch.index.query.QueryBuilders.queryStringQuery;

public class PatientHelper {

    public static RestHighLevelClient client;


    public static void main(String[] args) throws Exception {

        try (RestHighLevelClient clientTry = new RestHighLevelClient(
                RestClient.builder(
                        new HttpHost("localhost", 9200, "http")))) {

            client = clientTry;
            //
            List<String> listOfSpecimen = new ArrayList<>(Arrays.asList("SP00047", "SP00072", "SP00022"));
            List<String> list2 = new ArrayList<>(Arrays.asList("SP00011", "SP00061", "SP00036"));
            List<String> hpoTerms = getHpoTerms(getAPatientFromESFromID("PA00002"));
            hpoTerms.forEach(System.out::println);
            hpoTerms.clear();


            Map<String, Patient> donors = preparePedigree(listOfSpecimen);
            donors.forEach((k,v) -> {
                System.out.println("id="+k+"\n\t"+v);
//                System.out.println("\tspecimenId="+v.getSpecimenId());
//                System.out.println("\tFamilyid:"+v.getFamilyId());
//                System.out.println("\tisProband:"+v.isProban());
//                System.out.println("\tRelation:"+v.getRelation());
            });

            donors = preparePedigree(list2);
            donors.forEach((k,v) -> System.out.println("id=" + k + "\n\t" + v));

        }
    }


    static JSONObject getAPatientFromESFromID(String uid) {
        boolean requestSuccess = false;
        JSONObject patient = null;
        //BulkResponse bulkResponse;
        for (int i=0; i< 10; i++) {
            try {

                GetRequest getRequest = new GetRequest("patient", "patient", uid);

                GetResponse getResponse = client.get(getRequest, RequestOptions.DEFAULT);
                if (getResponse.isExists()) {

                    patient = new JSONObject(getResponse.getSourceAsString());

                }

                requestSuccess = true;
                break;
            } catch (Exception e) {
                System.err.println("****"+e);
            }
        }
        if (requestSuccess) {
//            System.out.println("*");
            return patient;
        } else {
            return null;
        }
    }

    static List<JSONObject> getPatientFromESWithQueryString(String query) {
        boolean requestSuccess = false;
        List<JSONObject> patients = new ArrayList<>();
        for (int i=0; i< 10; i++) {
            try {

                SearchRequest searchRequest = new SearchRequest("patient");
                SearchSourceBuilder searchSourceBuilder = new SearchSourceBuilder();

                QueryStringQueryBuilder queryStringQueryBuilder = queryStringQuery(query);
                searchSourceBuilder.query(queryStringQueryBuilder);
                searchRequest.source(searchSourceBuilder);

                SearchResponse response = client.search(searchRequest, RequestOptions.DEFAULT);

                Arrays.stream(response.getHits().getHits()).forEach((searchHit) -> {

                    JSONObject patient = new JSONObject(searchHit.getSourceAsString());
//                    System.out.println(searchHit.getSourceAsString());
                    patients.add(patient);

                });

                requestSuccess = true;
                break;
            } catch (Exception e) {
                System.err.println("****"+e);
            }
        }
        if (requestSuccess) {
//            System.out.println("*");
            return patients;
        } else {
            return null;
        }
    }


    static List<String> getHpoTerms(JSONObject patient) {


        // return a coma separated list of observed (POS) hpoTerms
        List<String> hpoTerms = new ArrayList<>();
        //JSONObject patient = getPatientFromES(patientId);
        //System.out.println("patient=\n"+patient.toString(0));

        if (patient == null) {
            return hpoTerms;
        }
        // Extract every phenotype code (for 1 patient) from his clinicalImpressions->observations->phenotypes if observed is POS

        JSONArray clinicalImpressions = (patient.isNull("clinicalImpressions") ? new JSONArray() : (JSONArray) patient.get("clinicalImpressions"));

        clinicalImpressions.forEach((ci) -> {

            JSONArray observations = ((JSONObject) ci).isNull("observations") ? new JSONArray() : (JSONArray) ((JSONObject) ci).get("observations");
            //System.out.println(((JSONObject) ci).toString());
            observations.forEach((obs) -> {
                if ( !((JSONObject) obs).isNull("phenotype") ) {
                    if ("POS".equalsIgnoreCase((String) ((JSONObject) obs).get("observed")) ) {
                        JSONArray phenotypes = (JSONArray) ((JSONObject) obs).get("phenotype");
                        phenotypes.forEach((phenotype) -> {
                            //System.out.println("\t" + ((JSONObject) phenotype).toString());
                            hpoTerms.add((String)  ((JSONObject) phenotype).get("code"));
                        });
                    }
                }
            });
        });

        return hpoTerms;
    }


    protected static String getId(JSONObject patient) {

        return (patient.isNull("id") ? "" : (String) patient.get("id"));
    }


    static JSONObject serviceRequest(JSONArray serviceRequests, String specimenId) {

        if (serviceRequests == null || serviceRequests.length() == 0) {
            return null;
        }

        for (int i=0; i<serviceRequests.length();i++) {

            JSONObject sr = (JSONObject) serviceRequests.get(i);
            JSONArray specimenArray = (JSONArray) sr.get("specimen");
            for (int j=0;j<specimenArray.length();j++) {
                JSONObject specimen = (JSONObject) specimenArray.get(j);
                if (specimenId.equalsIgnoreCase((String) specimen.get("id"))) {
                    return sr;
                }
            }
        }


        return null;
    }

    static String getPractitionerId(JSONArray practitioners, String requesterRoleId) {

        if (practitioners == null || practitioners.length() == 0) {
            return null;
        }

        for (int i=0; i<practitioners.length();i++) {

            JSONObject prac = (JSONObject) practitioners.get(i);

            String roleId = (String) prac.get("role_id");
            if (requesterRoleId.equalsIgnoreCase(roleId)) {
                return (String) prac.get("id");
            }

        }


        return null;
    }

    static String getOrgId(JSONArray practitioners, String requesterRoleId) {

        if (practitioners == null || practitioners.length() == 0) {
            return null;
        }

        for (int i=0; i<practitioners.length();i++) {

            JSONObject prac = (JSONObject) practitioners.get(i);

            String roleId = (String) prac.get("role_id");
            if (requesterRoleId.equalsIgnoreCase(roleId)) {
                return (String) prac.get("role_org_id");
            }

        }


        return null;
    }


    static Map<String, Patient> preparePedigree(List<String> listOfSpecimen) {
        Map<String, Patient> patientMap = new HashMap<>();

        listOfSpecimen.forEach((specimen) -> {
            Patient patient = new Patient(specimen);

            List<JSONObject> patients = getPatientFromESWithQueryString(specimen);
            JSONObject patientObj = patients.get(0);
            String studyId = (String) ((JSONObject) ((JSONArray) patientObj.get("studies")).get(0)).get("id");

            JSONObject serviceRequest = serviceRequest((JSONArray) patientObj.get("serviceRequests"), specimen);
            if (serviceRequest != null) {
                String code = (String) ((JSONObject) serviceRequest.get("code")).get("text");
                patient.setSequencingStrategy(code);
                String requesterRoleId = (String) serviceRequest.get("requester_id");
                // get practitionerId from requesterRoleId
                String practitionerId = getPractitionerId((JSONArray) patientObj.get("practitioners"), requesterRoleId);
                String orgId = getOrgId((JSONArray) patientObj.get("practitioners"), requesterRoleId);
                patient.setPractitionerId(practitionerId);
                patient.setOrgId(orgId);

            }

            patient.setStudyId(studyId);
            patient.setPatient(patientObj);
            String id = (String) patientObj.get("id");
            patient.setPatientId(id);
            patient.setFamilyId((String) patientObj.get("familyId"));
            boolean isProban = (Boolean) patientObj.get("isProband");
            patient.setProban(isProban);
            if (isProban) {
                patient.setRelation("Proban");
            }
            patient.setHposTerms(getHpoTerms(patientObj));
            patientMap.put(id, patient);
        });

        //hpoTerms.forEach((hpo) -> System.out.println(hpo));
        // relationship call 2
        if (listOfSpecimen.size() > 1 ) {
            patientMap.forEach((k,v) -> {
                if (v.isProban()) {
                    // fetch relationship of proban -- duo, trio, etc...
                    JSONArray link = (JSONArray) v.getPatient().get("link");
                    link.forEach((rel) -> {
                        String relationship = (String) ((JSONObject) rel).get("relationship");
                        String id = (String) ((JSONObject) rel).get("id");
                        patientMap.get(id).setRelation(relationship);
                    });
                }
            });
        }
        return patientMap;
    }


}
