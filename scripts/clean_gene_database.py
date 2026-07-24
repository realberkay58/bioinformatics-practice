import pandas as pd
import requests
import time

def has_clinical_data(gene_symbol):
    """MyVariant API'sine hızlıca bağlanıp genin ClinVar kaydı olup olmadığına bakar."""
    url = "https://myvariant.info/v1/query"
    query_str = f'(clinvar.gene.symbol:"{gene_symbol}" OR gene:"{gene_symbol}") AND _exists_:clinvar'
    
    try:
        response = requests.get(url, params={"q": query_str, "fields": "clinvar", "size": 1})
        if response.status_code == 200:
            data = response.json()
            return len(data.get("hits", [])) > 0
        return False
    except Exception:
        return False

def build_solid_database():
    print("🧬 Veritabanı temizliği başlıyor. Sadece klinik verisi olan genler seçilecek...")
    
    input_csv = "results/task05_differential_expression.csv" 
    output_csv = "results/task05_clinical_filtered.csv"
    
    try:
        df = pd.read_csv(input_csv)
    except FileNotFoundError:
        print(f"❌ {input_csv} bulunamadı! Terminalde doğru klasörde olduğundan emin ol.")
        return

    # Sütun adını otomatik bul (Gene, Symbol, ID vs.)
    gene_column = None
    for col in ["Gene", "gene", "Symbol", "Gene_Symbol", "ID", "id"]:
        if col in df.columns:
            gene_column = col
            break
            
    if not gene_column:
        print(f"❌ Gen isimlerinin olduğu sütun bulunamadı. Mevcut sütunlar: {df.columns.tolist()}")
        return
        
    print(f"✅ Gen sütunu bulundu: {gene_column}")
    
    valid_rows = []
    
    for index, row in df.iterrows():
        gene = str(row[gene_column]).strip()
        
        if has_clinical_data(gene):
            valid_rows.append(row)
            print(f"✅ {gene} eklendi. (Toplam Başarılı: {len(valid_rows)})")
            
        # 500 tane garantili gene ulaştığımızda dur
        if len(valid_rows) >= 500:
            print("🎉 500 adet klinik verisi olan gen başarıyla toplandı!")
            break
            
        time.sleep(0.1) # API'yi çökertmemek için kısa mola

    clean_df = pd.DataFrame(valid_rows)
    clean_df.to_csv(output_csv, index=False)
    print(f"📁 Temizlenmiş veritabanı şuraya kaydedildi: {output_csv}")

if __name__ == "__main__":
    build_solid_database()