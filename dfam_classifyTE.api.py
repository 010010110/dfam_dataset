import os
import re
import time
import requests

BASE_URL = "https://www.dfam.org/api"
OUTPUT_FASTA = "dfam_raw_target_elements.fasta"
LIMIT = 100


def is_model_target(classification, name=""):
    """Filtra se o elemento pertence ao conjunto de classes do modelo,

    sem alterar o cabeçalho original.
    """
    clean = re.sub(r"[;_/\-]", " ", classification).upper()
    tokens = set(clean.split())
    name_clean = re.sub(r"[;_/\-]", " ", name).upper()
    name_tokens = set(name_clean.split())

    # Descarta ruídos que não fazem parte do escopo (Crypton, SVA, snoRNA)
    if "CRYPTON" in clean or "SVA" in clean or "SNORNA" in clean:
        return False

    # 1.1 LTR: Copia, Gypsy, Bel-Pao
    if "COPIA" in tokens:
        return True
    if "GYPSY" in tokens and not any(
        v in clean for v in ["RETROVIRIDAE", "ORTHORETROVIRINAE", "MALR"]
    ):
        return True
    if "BEL" in tokens or "PAO" in tokens or "BEL PAO" in clean:
        return True

    # 1.2 DIRS
    if "DIRS" in tokens or "NGARO" in tokens:
        return True

    # 1.5 SINE: tRNA, 7SL, 5S ou SINE canônico
    if "SINE" in tokens or any(t in tokens for t in ["ALU", "MIR", "B2", "B4"]):
        return True

    # 1.4 LINE: L1, R2, RTE, Jockey, I (apenas autônomos legítimos)
    if "DEPENDENT" not in clean and "DERIVED" not in clean:
        line_subtypes = ["L1", "LINE1", "R2", "RTE", "JOCKEY"]
        if (
            any(t in tokens for t in line_subtypes)
            or any(t in name_tokens for t in line_subtypes)
            or "I ELEMENT" in clean
            or "LINE" in tokens
        ):
            return True

    # 2.1.1 TIRs / DNA Transposons
    is_class_ii = "CLASS II" in clean or "DNA TRANSPOSITION" in clean
    if is_class_ii:
        tir_targets = [
            "TC1",
            "MARINER",
            "TCMAR",
            "TC2",
            "TIGGER",
            "HAT",
            "MUTATOR",
            "MULE",
            "MUDR",
            "MERLIN",
            "TRANSIB",
            "PIGGYBAC",
            "HARBINGER",
            "PIF",
            "CACTA",
            "ENSPM",
            "CMC",
            "TIR",
            "TRANSPOSASE",
        ]
        if any(t in tokens for t in tir_targets) or any(
            t in name_tokens
            for t in ["MERLIN", "TRANSIB", "PIGGYBAC", "PB", "P"]
        ):
            return True
        if "P ELEMENT" in clean:
            return True

    return False


def get_target_families():
    print("Conectando à API do Dfam e selecionando famílias alvo...")
    headers = {"Accept": "application/json"}
    start = 1
    selected = []

    while True:
        params = {"start": start, "limit": LIMIT}
        try:
            r = requests.get(
                f"{BASE_URL}/families",
                params=params,
                headers=headers,
                timeout=30,
            )
            if r.status_code != 200:
                print(f"Erro {r.status_code} na página {start}")
                break

            data = r.json()
            results = data.get("results", [])
            if not results:
                break

            for item in results:
                acc = item.get("accession", "")
                # Considera apenas famílias curadas (DF...)
                if acc.startswith("DF"):
                    classification = item.get("classification") or ""
                    name = item.get("name") or ""
                    if is_model_target(classification, name):
                        selected.append(item)

            total_count = data.get("total_count", 0)
            print(
                f"  Processados {start + len(results) - 1}/{total_count} | Alvos selecionados: {len(selected)}"
            )

            if start + LIMIT > total_count or len(results) < LIMIT:
                break

            start += LIMIT
            time.sleep(0.04)

        except Exception as e:
            print(f"Exceção na consulta: {e}")
            break

    return selected


def download_fasta(accession):
    url = f"{BASE_URL}/families/{accession}/sequence"
    headers = {"Accept": "text/fasta"}
    params = {"format": "fasta"}
    try:
        r = requests.get(url, params=params, headers=headers, timeout=25)
        if r.status_code == 200:
            return r.text
    except Exception:
        pass
    return None


def main():
    targets = get_target_families()
    total = len(targets)
    print(f"\nTotal de famílias encontradas: {total}")

    if total == 0:
        print("Nenhum registro encontrado. Encerrando.")
        return

    with open(OUTPUT_FASTA, "w") as out:
        for idx, item in enumerate(targets, 1):
            acc = item["accession"]
            elem_name = item.get("name", "")
            classification = item.get("classification", "")

            seq_data = download_fasta(acc)
            if seq_data:
                lines = seq_data.strip().split("\n")
                if len(lines) > 1:
                    seq = "\n".join(lines[1:])
                    # Padrão nativo do Dfam: >Accession#Classificacao Nome
                    raw_header = f">{acc}#{classification} {elem_name}".strip()
                    out.write(f"{raw_header}\n{seq}\n")

            if idx % 50 == 0 or idx == total:
                print(f"  Progresso: {idx}/{total} sequências baixadas...")
            time.sleep(0.04)

    print(
        f"\nDownload concluído! Arquivo gerado com cabeçalhos brutos do Dfam em: '{OUTPUT_FASTA}'\n"
    )


if __name__ == "__main__":
    main()