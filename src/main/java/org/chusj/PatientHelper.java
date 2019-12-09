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

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
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
                //List<String> list2 = new ArrayList<>(Arrays.asList("SP00011", "SP00061", "SP00036"));
            List<String> hpoTerms = getHpoTerms(getAPatientFromESFromID("PA00002"));
            hpoTerms.forEach(System.out::println);
            hpoTerms.clear();

//            List<Pedigree> pedigrees = loadPedigree("pedigree.ped");
//            pedigrees.forEach((ped) -> System.out.println(ped));

            Map<String, Patient> donors = preparePedigree(listOfSpecimen);
            donors.forEach((k,v) -> {
                System.out.println("id="+k+"\n\t"+v);
//                System.out.println("\tspecimenId="+v.getSpecimenId());
//                System.out.println("\tFamilyid:"+v.getFamilyId());
//                System.out.println("\tisProband:"+v.isProban());
//                System.out.println("\tRelation:"+v.getRelation());
            });

            List<Pedigree> pedigrees = loadPedigree("pedigree.ped");
            donors = preparePedigreeFromPedAndFHIR(pedigrees);
            donors.forEach((k,v) -> System.out.println("id=" + k + "\n\t" + v));
            String pedigreePropsFile = "pedigree.properties";
            Properties pedigreeProps = VepHelper.getPropertiesFromFile(pedigreePropsFile);
            donors = preparePedigreeFromProps(pedigreeProps);
            System.out.println("--------");
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

    private static List<JSONObject> getPatientFromESWithQueryString(String query) {
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


    private static JSONObject serviceRequest(JSONArray serviceRequests, String specimenId) {

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

    private static String getPractitionerId(JSONArray practitioners, String requesterRoleId) {

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

    private static String getOrgId(JSONArray practitioners, String requesterRoleId) {

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
            boolean isInfected = (Boolean) patientObj.get("status");
            patient.setAffected(isInfected);
            patient.setGender((String) patientObj.get("gender"));
            patient.setProban(isProban);
            if (isProban) {
                patient.setRelation("Proban");
            }
            patient.setHposTerms(getHpoTerms(patientObj));
            patientMap.put(id, patient);
            patientMap.put(specimen, patient);
        });

        //hpoTerms.forEach((hpo) -> System.out.println(hpo));
        // relationship call 2
        if (listOfSpecimen.size() > 1 ) {
            patientMap.forEach((k,v) -> {
                if (v.isProban()) {
                    // fetch relationship of proban -- duo, trio, etc...

                    JSONArray link = (JSONArray) v.getPatient().get("link");

                    String[] linkToOthersPatientIds = new String[link.length()];
                    String[] linkToOthersSpecimenIds = new String[link.length()];
                    for (int i=0; i<link.length();i++) {
                        JSONObject rel = (JSONObject) link.get(i);
                        String relationship = (String) rel.get("relationship");
                        String id = (String) rel.get("id");
                        linkToOthersPatientIds[i] = id;
                        System.out.println("V:"+v.getPatientId()+"id:"+id);
                        linkToOthersSpecimenIds[i] = patientMap.get(id).getSpecimenId();
                        patientMap.get(id).setRelation(relationship);

                    }
                    v.setLinkToOthersPatientId(linkToOthersPatientIds);
                    v.setLinkToOthersSpecimenId(linkToOthersSpecimenIds);
                }
            });
        }
        return patientMap;
    }

    static Map<String, Patient> preparePedigreeFromPedAndFHIR(List<Pedigree> pedigrees) {
//        List<Pedigree> pedigrees = loadPedigree(filename);

        List<String> listOfSpecimen = new ArrayList<>();
        pedigrees.forEach((ped) -> listOfSpecimen.add(ped.getId()));

        //listOfSpecimen.forEach(System.out::println);
        return preparePedigree(listOfSpecimen);

    }

    /*
    #into pedigree, put specimen IDs in comma separated list

    linkQty=2,0,0
    linkPatient=PA00027,PA00052,.,.,.,.
    linkSpecimen=14141,14142,.,.,.,.
    patientId=PA00002,PA00027,PA00052
    familyId=FA0002,FA0002,FA0002
    studyId=ET00011,ET00011,ET00011
    organizationId=OR00202,OR00202,OR00202

    annotationTool=VEP 97
    HpoQtyTerms=2,0,0
    hpoTermsPos=HP:0005280,HP:0001773

    */
    public static Map<String, Patient> preparePedigreeFromProps(Properties pedigreeProps) {
        Map<String, Patient> patientMap = new HashMap<>();

        String[] specimens = pedigreeProps.getProperty("specimen").split(",");
        String[] familyIds = pedigreeProps.getProperty("familyId").split(",");
        String[] genders = pedigreeProps.getProperty("gender").split(",");
        String[] patientIds = pedigreeProps.getProperty("patientId").split(",");
        String[] relationship = pedigreeProps.getProperty("relation").split(",");
        String[] studyIds = pedigreeProps.getProperty("studyId").split(",");
        String sequencingStrategy = pedigreeProps.getProperty("sequencingStrategy");
        String[] organizationId = pedigreeProps.getProperty("organizationId").split(",");
        String[] practitionerId = pedigreeProps.getProperty("practitionerId").split(",");
        String[] laboNames = pedigreeProps.getProperty("laboName").split(",");
        String[] isInfected = pedigreeProps.getProperty("isAffected").split(",");
        String[] isProban = pedigreeProps.getProperty("isProban").split(",");
        String[] hpoTermsPos = pedigreeProps.getProperty("hpoTermsPos").split(",");
        String[] hpoQtyTerms = pedigreeProps.getProperty("HpoQtyTerms").split(",");
        String[] linkQty = pedigreeProps.getProperty("linkQty").split(",");
        String[] linkPatient = pedigreeProps.getProperty("linkPatient").split(",");
        String[] linkSpecimen = pedigreeProps.getProperty("linkSpecimen").split(",");

        int hpoArrayPosition=0;
        int pos =0;
        for (int i=0; i<specimens.length;i++) {
            Patient patient = new Patient(specimens[i]);
            patient.setGender(genders[i]);
            patient.setAffected(Boolean.parseBoolean(isInfected[i]));
            patient.setRelation(relationship[i]);
            patient.setFamilyId(familyIds[i]);

            int hpoTermsQty = Integer.parseInt(hpoQtyTerms[i]);
            List<String> hpoList = new ArrayList<>();
            for (int j=0;j<hpoTermsQty;j++) {
                hpoList.add(hpoTermsPos[hpoArrayPosition++]);
            }
            patient.setHposTerms(hpoList);
            patient.setOrgId(organizationId[i]);
            patient.setPatientId(patientIds[i]);
            patient.setPractitionerId(practitionerId[i]);
            patient.setStudyId(studyIds[i]);
            patient.setLabName(laboNames[i]);
            patient.setProban(Boolean.parseBoolean(isProban[i]));
            patient.setSequencingStrategy(sequencingStrategy);

            int linkNumber = Integer.parseInt(linkQty[i]);
            String[] linkPatientTo = new String[linkNumber];
            String[] linkSpecimenTo = new String[linkNumber];
            for (int j=0;j<linkNumber;j++) {
                linkPatientTo[j] = linkPatient[pos];
                linkSpecimenTo[j] = linkSpecimen[pos++];
            }
            patient.setLinkToOthersPatientId(linkPatientTo);
            patient.setLinkToOthersSpecimenId(linkSpecimenTo);

            patientMap.put(specimens[i],patient);
            patientMap.put(patientIds[i],patient);
        }

        return patientMap;
    }


    public static List<Pedigree> loadPedigree(String filename) {
        List<Pedigree> pedigrees = new ArrayList<>();

        try (BufferedReader buf = new BufferedReader(new FileReader(filename))) {

            String fetchedLine;
            while (true) {
                fetchedLine = buf.readLine();
                if (fetchedLine == null || fetchedLine.length() == 0) {
                    break;
                } else {
                    Pedigree ped = new Pedigree();
                    String[] line = fetchedLine.split("\t");
                    //System.out.println(fetchedLine + "-"+line.length);
                    ped.setFamilyId(line[0]);
                    ped.setId(line[1]);
                    ped.setPaternalId(line[2]);
                    ped.setMaternalId(line[3]);
                    ped.setSex(line[4]);
                    ped.setPhenotype(line[5]);
                    pedigrees.add(ped);

                }
            }

        } catch (IOException e) {
            e.printStackTrace();
        }

        return pedigrees;
    }


}
