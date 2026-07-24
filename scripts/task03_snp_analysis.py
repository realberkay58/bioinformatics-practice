import requests

def fetch_clinvar_myvariant(gene_symbol):
    print(f"📡 Fetching clinical variants via Ensembl REST API for: {gene_symbol}")
    
    server = "https://rest.ensembl.org"
    # Ensembl overlap API: Gen sembolüyle fiziksel olarak çakışan tüm varyantları getirir
    ext = f"/overlap/symbol/human/{gene_symbol}?feature=variation"
    
    # Ensembl sunucusuna veriyi JSON formatında istediğimizi belirtiyoruz
    headers = { "Content-Type" : "application/json" }
    
    try:
        response = requests.get(server + ext, headers=headers)
        
        if response.status_code == 200:
            data = response.json()
            variants_list = []
            
            for variant in data:
                # Sadece klinik önemi (clinical_significance) laboratuvarca raporlanmışları filtrele
                clin_sig = variant.get('clinical_significance', [])
                
                if clin_sig:
                    # Ensembl bu veriyi dizi (array) döner. String'e çevirip baş harflerini büyütüyoruz (Pathogenic, vs.)
                    significance = ", ".join(clin_sig).title()
                    
                    # Ensembl ID olarak genelde evrensel rsID (dbSNP) numarasını döner
                    variant_id = variant.get('id', 'Unknown ID')
                    
                    # Mutasyon lokasyonu: Ensembl HGVS yerine Kromozom ve Allel koordinatı verir (Örn: Chr17:7577538 (C/T))
                    alleles = variant.get('alleles', 'Unknown')
                    seq_region = variant.get('seq_region_name', '')
                    start = variant.get('start', '')
                    coding = f"Chr{seq_region}:{start} ({alleles})"
                    
                    # Ensembl 'overlap' endpoint'i hastalık ismini doğrudan dönmez
                    condition = "Check NCBI/dbSNP for exact condition"
                    
                    variants_list.append({
                        "Variant": coding,
                        "ClinVar ID": variant_id,
                        "Significance": significance,
                        "Condition": condition
                    })
                    
                    # Sayfa hızını korumak için en önemli 15 varyantta kesiyoruz
                    if len(variants_list) >= 15:
                        break
                        
            return variants_list
        else:
            return []
            
    except Exception as e:
        print(f"❌ Ensembl API Error: {e}")
        return []