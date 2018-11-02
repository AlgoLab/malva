import sys
import csv, statistics

def vcf_analysis(cvs_content):
    # Analysis on samples (original vs corrected)
    res_dict = {}
    for row in cvs_content:
        row_key = "_".join([row[0],row[1],row[2]])
        p,r = float(row[-4]), float(row[-3])
        res_dict[row_key] = res_dict[row_key] + [(row[3],p,r)] if row_key in res_dict else [(row[3],p,r)]

    dps = [] # differences between precision values
    drs = [] # differences between recall values
    for (k,[(t1,p1,r1), (t2,p2,r2)]) in res_dict.items():
        if t1 == "full":
            dps += [round(p1-p2,3)]
            drs += [round(r1-r2,3)]
        else:
            dps += [round(p2-p1,3)]
            drs += [round(r2-r1,3)]
        #print(dps[-1], drs[-1])
    #print("---")
    print(round(statistics.mean(dps),3), round(statistics.mean(drs),3))

def samples_analysis(cvs_content):
    # Analysis on samples (original vs corrected)
    res_dict = {}
    for row in cvs_content:
        row_key = "_".join([row[0],row[1],row[3]])
        p,r = float(row[-4]), float(row[-3])
        res_dict[row_key] = res_dict[row_key] + [(row[2],p,r)] if row_key in res_dict else [(row[2],p,r)]

    dps = [] # differences between precision values
    drs = [] # differences between recall values
    for (k,[(t1,p1,r1), (t2,p2,r2)]) in res_dict.items():
        if t1 == "corrected":
            dps += [round(p1-p2,3)]
            drs += [round(r1-r2,3)]
        else:
            dps += [round(p2-p1,3)]
            drs += [round(r2-r1,3)]
        #print(dps[-1], drs[-1])
    #print("---")
    print(round(statistics.mean(dps),3), round(statistics.mean(drs),3))

def main():
    csv_path = sys.argv[1]
    csv_file = open(csv_path, 'r', newline='')
    csv_content = csv.reader(csv_file, delimiter=',', quotechar='|')
    next(csv_content, None) # skip the headers

    samples_analysis(csv_content)

    # csv_file.seek(0)
    # next(csv_content, None) # skip the headers
    # vcf_analysis(csv_content)

if __name__ == '__main__':
    main()
