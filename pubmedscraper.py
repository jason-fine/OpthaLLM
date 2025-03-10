import os
from Bio import Entrez

# Set your email (required by NCBI API)
Entrez.email = "keuriglover949@gmail.com"

# File to store retrieved PubMed IDs
ID_FILE = "retrieved_ids.txt"

def load_retrieved_ids():
    """Load previously retrieved PubMed IDs from a file."""
    if os.path.exists(ID_FILE):
        with open(ID_FILE, "r") as f:
            return set(f.read().splitlines())
    return set()

def save_retrieved_ids(pubmed_ids):
    """Save newly retrieved PubMed IDs to a file."""
    with open(ID_FILE, "a") as f:
        for pubmed_id in pubmed_ids:
            f.write(pubmed_id + "\n")

def search_pubmed(query, max_results=10):
    """Search PubMed for articles matching the query, only fetching new ones."""
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="pub_date")
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

def fetch_article_details(pubmed_ids):
    """Fetch article details using PubMed IDs."""
    if not pubmed_ids:
        return []
    
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
    max_results = 5

    retrieved_ids = load_retrieved_ids()
    new_pubmed_ids = search_pubmed(query, max_results=max_results)

    # Filter out already retrieved articles
    new_pubmed_ids = [pid for pid in new_pubmed_ids if pid not in retrieved_ids]

    if new_pubmed_ids:
        articles = fetch_article_details(new_pubmed_ids)
        for article in articles:
            print(f"Title: {article['title']}")
            print(f"Authors: {', '.join(article['authors'])}")
            print(f"Abstract: {article['abstract']}\n{'-'*80}")

        # Save newly retrieved IDs
        save_retrieved_ids(new_pubmed_ids)
    else:
        print("No new articles found today.")
