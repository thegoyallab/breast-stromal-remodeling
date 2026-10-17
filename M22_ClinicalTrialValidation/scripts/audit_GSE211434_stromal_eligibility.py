import gzip,csv,collections
p='/mnt/stk35data/BreastCancer_Census_v1/StageM22_ClinicalTrialStromalValidation_v0_1/raw/GSE211434_cell_types.txt.gz'
rows=[]
with gzip.open(p,'rt',newline='') as f:
    rd=csv.DictReader(f,delimiter='\t')
    rows.extend(rd)
print('cells',len(rows))
patients=sorted(set(r['Patient'] for r in rows))
print('patients',len(patients))
print('days',dict(collections.Counter(r['Day'] for r in rows)))
print('celltypes',collections.Counter(r['CellType'] for r in rows).most_common())
targets=['Fibroblasts','Pericytes','Endothelial cells']
cnt=collections.Counter((r['Patient'],int(r['Day']),r['CellType']) for r in rows)
meta={}
for r in rows:
    meta[r['Patient']]={'responseYN':r['ResponseYN'],'response':r['Response'],'treatment':r['Treatment'],'arm':r['TreatmentArm']}
for t in targets:
    print('LINEAGE',t)
    for thr in [20,10]:
        ev={day:sum(cnt[(pt,day,t)]>=thr for pt in patients) for day in [0,14,180]}
        p14=sum(cnt[(pt,0,t)]>=thr and cnt[(pt,14,t)]>=thr for pt in patients)
        p180=sum(cnt[(pt,0,t)]>=thr and cnt[(pt,180,t)]>=thr for pt in patients)
        triple=sum(all(cnt[(pt,d,t)]>=thr for d in [0,14,180]) for pt in patients)
        print('thr',thr,'evaluables',ev,'paired0_14',p14,'paired0_180',p180,'triple',triple)
        for day in [14,180]:
            cc=collections.Counter()
            for pt in patients:
                if cnt[(pt,0,t)]>=thr and cnt[(pt,day,t)]>=thr:
                    cc[meta[pt]['responseYN']]+=1
            print('response_by_pair','day',day,'thr',thr,dict(cc))
out='/mnt/stk35data/BreastCancer_Census_v1/StageM22_ClinicalTrialStromalValidation_v0_1/derived/GSE211434_patient_day_lineage_cell_counts.tsv'
celltypes=sorted(set(r['CellType'] for r in rows))
with open(out,'w',newline='') as f:
    w=csv.writer(f,delimiter='\t')
    w.writerow(['Patient','Day','Treatment','TreatmentArm','Response','ResponseYN','CellType','n_cells'])
    for pt in patients:
        for d in [0,14,180]:
            for t in celltypes:
                n=cnt[(pt,d,t)]
                if n:
                    m=meta[pt]
                    w.writerow([pt,d,m['treatment'],m['arm'],m['response'],m['responseYN'],t,n])
print('wrote',out)
