import pandas as pd
import os

def process_reports(directory, classification_code):
    data = {}
    report_files = [f for f in os.listdir(directory) if f.endswith(".report")]

    for report_file in report_files:
        file_path = os.path.join(directory, report_file)
        with open(file_path, "r") as f:
            lines = f.readlines()
        
        taxon_data = {}
        for line in lines:
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            
            percentage, total_fragments, direct_fragments, code, taxon_id, taxon_name = parts[:6]
            
            # Se o número de fragmentos diretos for 0, usa o total de fragmentos cobertos pelo clado
            direct_fragments = int(direct_fragments) if int(direct_fragments) > 0 else int(total_fragments)

            if code == classification_code:
                taxon_data[taxon_id] = direct_fragments

        data[report_file] = taxon_data

        # Criando um dataframe específico para este arquivo
        df_single = pd.DataFrame(list(taxon_data.items()), columns=["Taxon ID", "Fragmentos"])
        txt_filename_single = f"{report_file}_{classification_code}_dataframe.txt"
        df_single.to_csv(os.path.join(directory, txt_filename_single), sep="\t", index=False)

    # Criando o DataFrame consolidado
    df_total = pd.DataFrame.from_dict(data, orient="index").transpose()
    df_total.fillna(0, inplace=True)

    # Salvando o dataframe consolidado como .txt
    txt_filename_total = f"resultado_total_{classification_code}_dataframe.txt"
    df_total.to_csv(os.path.join(directory, txt_filename_total), sep="\t", index=True)

    return df_total

# Uso:
directory = "~/Documents/kraken_carine"  # Ajuste para o seu diretório
classification_code = "F"  # Código de classificação desejado
df_final = process_reports(os.path.expanduser(directory), classification_code)

print("Arquivos salvos com sucesso!")
