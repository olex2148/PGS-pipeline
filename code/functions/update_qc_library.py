#!/usr/bin/env python

"""Update the QC library



"""
 
import csv, sys
from pathlib import Path

from openpyxl import load_workbook

qc_file, library_folder, library_name, run_date = sys.argv[1:]
library_dir = Path(library_folder)
#Lock the file while is being filled, we are going to write it while each job is running
library_file = max(library_dir.glob(f"????_??_??_{library_name}"), default=library_dir / library_name)
output_file = str(library_dir / f"{run_date}_{library_name}")

with open(qc_file, newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    source_columns = reader.fieldnames
    #Parameters we keep
    value_columns = [
        "h2_init", "h2_init_se", "h2_init_z", "h2_auto", "h2_2.5", "h2_97.5",
        "n_chains", "shrinkage", "p_auto", "ldsc_ratio"
        ]
        
    selected_columns = [
      column for column in source_columns
      if column in value_columns
      or column.endswith("_qc")
      ]
    if "ID" not in source_columns:
        raise ValueError("No ID collumn in QC results")
          
    #Convert to numbers for Excel
    def excel_value(value):
        if value in (None,"","NA"):
            return None
        try:
            number = float(value)
            return number if number not in (float("inf"), float("-inf")) else value
        except ValueError:
            return value
        
    #Explain simply on cases of WARN and FAIL
    def make_qc_comments(row): #Nothing in PASS
        if row.get("overall_qc") == "PASS":
            return None        
          
        #If WARN/FAIL add commentary
        comments = []
        def value(parameter):
            raw = row.get(parameter)
            try:
                return f"{float(raw):.3g}"
            except (TypeError, ValueError):
                return str(raw)
       
        for column in selected_columns:
            if column == "overall_qc" or not column.endswith("_qc"):
                continue
               
            flag = row.get(column, "")
            if not flag.startswith(("WARN","FAIL")):
                continue
               
            parameter = column[:-3]
               
            parameter ={
                "h2_power": "h2_init_z",
                "chains": "n_chains",
            }.get(parameter, parameter)
               
            if "missing" in flag.lower(): comment = f"{parameter} missing"
            elif column == "h2_init_qc": comment = ("Implausible LDSC heritability estimate " f"(h2 = {value("h2_init")}).")
            elif column == "h2_power_qc": comment = (f"Low LDSC heritability Z ({value("h2_init_z")}).")
            elif column == "h2_auto_qc": comment = ("Implausible LDpred2 heritability estimate " f"(h2 = {value("h2_auto")})")
            elif column == "chains_qc": comment = (f"Only {value("n_chains")}/50 chains passed")
            elif column == "shrinkage_qc": comment = (f"Strong shrinkage ({value("shrinkage")})")
            elif column == "p_auto_qc": comment = (f"Polygenicity outside (0,1]({value("p_auto")})")
            elif column == "ldsc_ratio_qc": 
                comment =(f"ldsc ratio ({value("ldsc_ratio")})")
                if flag.startswith("FAIL"):
                    comment += ( " suggests possible counfunding or population structure")
            else: comment = f"{parameter}: {flag}"
            comments.append(comment)
        return " ".join(comments) or None
           
        
    new_rows =[]
        
    for row in reader:
        pgs_code = row["ID"].strip()
        result ={"PGS_code": pgs_code}
        for column in selected_columns:
            result[column] = excel_value(row[column])
                
        result["Date"] = run_date
        result["QC_comments"] = make_qc_comments(row)
        new_rows.append(result)
       

#Adding to the library
workbook = load_workbook(library_file)
sheet = workbook["06_PGS_QC"]
#Read previous rows
existing = list(sheet.iter_rows(values_only=True))
old_columns = (list(existing[0])
    if existing and any(value is not None for value in existing[0])
    else []
)
old_rows = [dict(zip(old_columns, row))
    for row in existing[1:]
    if any(value is not None for value in row)
]

#Adding the new columns to the existing ones
columns = ["PGS_code"]
updated_traits ={row["PGS_code"] for row in new_rows}
retained_rows =[ row for row in old_rows if str(row.get("PGS_code", "")).strip() not in updated_traits]
for column in old_columns + selected_columns:
    if column and column not in columns and column not in ("Date", "overall_qc", "QC_comments"): columns.append(column)
columns.append("Date")
    
if "overall_qc" in old_columns or "overall_qc" in selected_columns: columns.append("overall_qc")
columns.append("QC_comments")
#Rewrite the QC sheet
sheet.delete_rows(1, sheet.max_row)
sheet.append(columns)
    
for row in retained_rows + new_rows: sheet.append([row.get(column) for column in columns])
        
sheet.freeze_panes = "A2"
sheet.auto_filter.ref = sheet.dimensions
    
#Save to a temporary file in case it breaks suddenly
output_path = Path(output_file)
temporary_path = output_path.with_name(f".{output_path.stem}.tmp.xlsx")
  
try:
    workbook.save(temporary_path)
    temporary_path.replace(output_path)
finally:
    workbook.close()
    temporary_path.unlink(missing_ok=True)
