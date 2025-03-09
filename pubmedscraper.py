import requests
from bs4 import BeautifulSoup
from Bio import Entrez

# Set your email (required by NCBI API)
Entrez.email = "keuriglover949@gmail.com"

def search_pubmed(query, max_results=10):
    """Search PubMed for articles matching the query."""
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

def fetch_article_details(pubmed_ids):
    """Fetch article details using PubMed IDs."""
    handle = Entrez.efetch(db="pubmed", id=",".join(pubmed_ids), retmode="xml")
    articles = Entrez.read(handle)
    handle.close()
    
    results = []
    for article in articles["PubmedArticle"]:
        title = article["MedlineCitation"]["Article"]["ArticleTitle"]
        authors = [
            author["LastName"] + " " + author["ForeName"]
            for author in article["MedlineCitation"]["Article"].get("AuthorList", [])
        ]
        abstract = article["MedlineCitation"]["Article"].get("Abstract", {}).get("AbstractText", ["No abstract"])[0]
        results.append({"title": title, "authors": authors, "abstract": abstract})
    
    return results

if __name__ == "__main__":
    query = "ophthalmology"
    pubmed_ids = search_pubmed(query, max_results=5)  # Get 5 article IDs
    articles = fetch_article_details(pubmed_ids)  # Fetch article details
    
    for article in articles:
        print(f"Title: {article['title']}")
        print(f"Authors: {', '.join(article['authors'])}")
        print(f"Abstract: {article['abstract']}\n{'-'*80}")
