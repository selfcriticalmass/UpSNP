import urllib.request
import urllib.parse
import json
import time
import re
import sys

EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
VARIATION_URL = "https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/"

def chunk_list(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def normalize(raw):
    return re.sub(r'^rs', '', raw.strip(), flags=re.IGNORECASE)

def fetch_json(url, params=None, data=None):
    headers = {"User-Agent": "snptracker_run/0.1.0"}
    if data:
        data = urllib.parse.urlencode(data).encode()
    if params:
        url = url + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, data=data, headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            return json.loads(resp.read().decode())
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return None
        print(f"  HTTP error {e.code}: {e.reason}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"  Network error: {e}", file=sys.stderr)
        return None

def resolve_via_entrez(rsid_list):
    results = {}
    for chunk in chunk_list(rsid_list, 200):
        nums = [normalize(r) for r in chunk]
        print(f"  Querying {len(nums)} rsIDs via E-utilities...", file=sys.stderr)
        data = fetch_json(EUTILS_URL, params={"db": "snp", "id": ",".join(nums), "retmode": "json"})
        if data and "result" in data:
            for num in nums:
                qk = f"rs{num}"
                record = data["result"].get(num)
                if record:
                    can_snp = record.get("snp_id")
                    results[qk] = {
                        "latest_id": f"rs{can_snp}" if can_snp else "Not_Found_Or_Deleted",
                        "coords": record.get("chrpos", "Unknown"),
                    }
                else:
                    results[qk] = {"latest_id": "Not_Found_Or_Deleted", "coords": "Unknown"}
        else:
            for num in nums:
                results[f"rs{num}"] = {"latest_id": "Not_Found_Or_Deleted", "coords": "Unknown"}
        time.sleep(0.3)
    return results

def _extract_spdi_coords(data):
    placements = data.get("primary_snapshot_data", {}).get("placements_with_allele", [])
    for p in placements:
        seq = p.get("seq_id", "")
        if seq.startswith("NC_"):
            alleles = p.get("alleles", [])
            if alleles:
                spdi = alleles[0].get("allele", {}).get("spdi", {})
                pos = spdi.get("position")
                if pos is not None:
                    return f"{seq}:{pos}"
    return "Unknown"

def resolve_via_variation_api(rsid_list):
    results = {}
    for raw in rsid_list:
        num = normalize(raw)
        qk = f"rs{num}"
        data = fetch_json(f"{VARIATION_URL}{num}")
        if data is None:
            results[qk] = {"latest_id": "Not_Found_Or_Deleted", "coords": "Unknown"}
        else:
            merged = data.get("merged_snapshot_data", {})
            if merged.get("merged_into"):
                latest_num = merged["merged_into"][0]
                latest = f"rs{latest_num}"
                current_data = fetch_json(f"{VARIATION_URL}{latest_num}")
                coords = _extract_spdi_coords(current_data) if current_data else "Unknown"
            else:
                latest = qk
                coords = _extract_spdi_coords(data)

            results[qk] = {"latest_id": latest, "coords": coords}
        time.sleep(0.5)
    return results

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Map legacy rsIDs to current canonical rsIDs with coordinates via NCBI APIs")
    parser.add_argument("-i", "--input", default="historical_rsids.txt",
                        help="Input file, one rsID per line")
    parser.add_argument("-o", "--output", default="updated_rsids_manifest.tsv",
                        help="Output TSV file")
    parser.add_argument("--method", choices=["entrez", "variation"], default="entrez",
                        help="API: 'entrez' (E-utilities, batch) or 'variation' (per-ID)")
    parser.add_argument("--assembly", default=None,
                        help="Filter to a specific RefSeq assembly accession prefix (e.g. NC_000011 for chr11)")
    args = parser.parse_args()

    with open(args.input) as f:
        legacy_ids = [line.strip() for line in f if line.strip()]

    if not legacy_ids:
        print("No rsIDs found.", file=sys.stderr)
        sys.exit(1)

    print(f"Loaded {len(legacy_ids)} rsIDs from {args.input}", file=sys.stderr)
    print(f"Method: {args.method}", file=sys.stderr)

    resolver = resolve_via_entrez if args.method == "entrez" else resolve_via_variation_api
    updated_map = resolver(legacy_ids)

    with open(args.output, "w") as out:
        out.write("Legacy_rsID\tLatest_rsID\tCoordinates\tStatus\n")
        for old_id in legacy_ids:
            norm = normalize(old_id)
            qk = f"rs{norm}"
            entry = updated_map.get(qk, {"latest_id": "Not_Found_Or_Deleted", "coords": "Unknown"})
            new_id = entry["latest_id"]
            coords = entry.get("coords", "Unknown")
            if new_id == qk:
                status = "Unchanged"
            elif new_id == "Not_Found_Or_Deleted":
                status = "Withdrawn"
            else:
                status = "Merged"
            out.write(f"{old_id}\t{new_id}\t{coords}\t{status}\n")

    print(f"Done -> {args.output}", file=sys.stderr)
