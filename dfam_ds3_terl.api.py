import time
import requests

BASE_URL = "https://www.dfam.org/api"

# Termos-alvo para casar com a classificação/linhagem
TARGET_TERMS = [
    "Copia",
    "Gypsy",
    "Bel-Pao",
    "Pao",
    "ERV",
    "L1",
    "SINE",
    "Tc1-Mariner",
    "TcMar",
    "hAT",
    "Mutator",
    "MULE",
    "PIF-Harbinger",
    "Harbinger",
    "CACTA",
    "EnSpm",
]

# TARGET_TERMS = [
#     "Copia",
#     "Gypsy",
#     "Bel-Pao",
#     "Pao",
#     "ERV",
#     "L1",
#     "Tc1-Mariner",
#     "TcMar",
#     "hAT"
# ]

output_file = "dfam_benchmark_ds1-terl.fasta"


def get_curated_families():
    print("Conectando à API do Dfam...")
    families = []
    start = 1
    limit = 100

    headers = {"Accept": "application/json"}

    while True:
        url = f"{BASE_URL}/families"
        params = {"start": start, "limit": limit}

        response = requests.get(url, params=params, headers=headers, timeout=30)
        if response.status_code != 200:
            print(f"Erro {response.status_code}: {response.text}")
            break

        data = response.json()
        results = data.get("results", [])
        if not results:
            break

        for item in results:
            acc = item.get("accession", "")
            # DF = Curated / DR = Uncurated (repetições cruas de genoma)
            if acc.startswith("DF"):
                families.append(item)

        total_count = data.get("total_count", 0)
        print(
            f"  Progresso: {len(families)} curadas encontradas (processando {start} de {total_count})..."
        )

        if start + limit > total_count or len(results) < limit:
            break

        start += limit
        time.sleep(0.05)

    return families


def download_fasta(accession):
    url = f"{BASE_URL}/families/{accession}/sequence"
    # O endpoint de sequência aceita format=fasta ou Header Accept
    headers = {"Accept": "text/fasta"}
    params = {"format": "fasta"}
    try:
        r = requests.get(url, params=params, headers=headers, timeout=20)
        if r.status_code == 200:
            return r.text
    except Exception as e:
        print(f"Erro em {accession}: {e}")
    return None


# 1. Obter todas as famílias curadas
all_curated = get_curated_families()
print(f"\nTotal de famílias curadas recuperadas: {len(all_curated)}")

# 2. Filtrar pelas linhagens desejadas
selected = []
for f in all_curated:
    classification = f.get("classification") or ""
    if any(term.lower() in classification.lower() for term in TARGET_TERMS):
        selected.append(f)

print(f"Total selecionado para download: {len(selected)}")

# 3. Baixar os FASTAs
with open(output_file, "w") as out:
    for idx, f in enumerate(selected, 1):
        acc = f.get("accession")
        name = f.get("name", "")
        classification = f.get("classification", "")

        print(f"[{idx}/{len(selected)}] Baixando {acc} ({classification})...")
        fasta_content = download_fasta(acc)

        if fasta_content:
            lines = fasta_content.strip().split("\n")
            if len(lines) > 1:
                seq = "\n".join(lines[1:])
                header = f">{acc}#{classification} {name}".strip()
                out.write(f"{header}\n{seq}\n")

        time.sleep(0.05)

print(f"\nArquivo salvo com sucesso em '{output_file}'.")