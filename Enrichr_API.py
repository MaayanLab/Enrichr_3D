
import json
import requests

## Analyze gene list
def enrichr_submitQuery(geneList, description):
    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/addList'
    genes_str = '\n'.join(geneList)
    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }
    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')
    ID_info = json.loads(response.text)
    # print(ID_info)
    return ID_info

## View added gene list
def enrichr_viewGeneList(ID_info):
    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/view?userListId=%s'
    user_list_id = ID_info['userListId']#363320
    response = requests.get(ENRICHR_URL % user_list_id)
    if not response.ok:
        raise Exception('Error getting gene list')
    gene_list = json.loads(response.text)
    # print(gene_list)
    return gene_list

## Get enrichment results
def enrichr_enrichmentResults(ID_info, gene_set_library = 'KEGG_2015'):
    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'
    user_list_id = ID_info['userListId']#363320
    response = requests.get(
        ENRICHR_URL + query_string % (user_list_id, gene_set_library)
     )
    if not response.ok:
        raise Exception('Error fetching enrichment results')
    results = json.loads(response.text)
    # print(results)
    return results

## Find terms that contain a given gene
def enrichr_findGeneTerms(gene):
    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/genemap'
    query_string = '?json=true&setup=true&gene=%s'
    gene = gene#'AKT1'
    response = requests.get(ENRICHR_URL + query_string % gene)
    if not response.ok:
        raise Exception('Error searching for terms')
    gene_search = json.loads(response.text)
    print(gene_search)
    return gene_search

## Download file of enrichment results
def enrichr_downloadResults(ID_info):
    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/export'
    query_string = '?userListId=%s&filename=%s&backgroundType=%s'
    user_list_id = ID_info['userListId']#363320
    filename = 'example_enrichment'
    gene_set_library = 'KEGG_2015'
    url = ENRICHR_URL + query_string % (user_list_id, filename, gene_set_library)
    response = requests.get(url, stream=True)
    with open(filename + '.txt', 'wb') as f:
        for chunk in response.iter_content(chunk_size=1024):
            if chunk:
                f.write(chunk)

# Run Enrichr for each cluster
def run_Enrichr(clusterDict):
    enrichrResults={}
    for clustID in clusterDict:
        # print("Processing : "+ str(clustID))
        cluster = clusterDict[clustID]
        geneList = cluster['predictedKinases']+cluster['targetKinases']
        ID_info = enrichr_submitQuery(geneList=geneList, description=str(clustID) )
        results='Fail'
        while results=='Fail':
            try:
                results = enrichr_enrichmentResults(ID_info=ID_info, gene_set_library = 'KEGG_2015')
                enrichrResults[clustID] = {'results':results,'geneList':geneList}
            except:
                results = 'Fail'
    # Get the top terms for each cluster
    names=['Rank', 'Term name', 'P-value', 'Z-score', 'Combined score', 'Overlapping genes', 'Adjusted p-value', 'Old p-value',
           'Old adjusted p-value']
    enrichrDF = pd.DataFrame([enrichrResults[x]['results']['KEGG_2015'][0] for x in enrichrResults], columns=names, index=None)
    enrichrDF['geneList length'] = [len(enrichrResults[x]['geneList']) for x in enrichrResults]
    enrichrDF['Library'] = 'KEGG_2015'
    enrichrDF['clusterID'] = enrichrResults.keys()
    return enrichrDF


# # Raw Ranks
# zRanks_enrichrDF = run_Enrichr(zRanks_clustDict)
# zRanks_enrichrDF.sort_values(by='geneList length',ascending=False, inplace=True)
#
# # Correlation Matrix
# corr_enrichrDF = run_Enrichr(corr_clustDict)
#
#
#
# corr_enrichrSig = corr_enrichrDF[corr_enrichrDF['Adjusted p-value']<=0.05]
# print(str(round(len(corr_enrichrSig)/len(corr_enrichrDF)*100, 2)) + \
#       "% of the cluster modules are significantly enriched for ontological terms.")
# corr_enrichrSig.head()