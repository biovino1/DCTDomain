import os

def main():

    linkage = 'ward'
    metric = 'euclidean'
    linkage_data = 'data/linkage_data.pkl'
    thresholds = [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000]
    criterion = 'distance'

    for t in thresholds:
        os.system(f'python cluster_dct.py -m {linkage} -p {metric} -l {linkage_data} -t {t} -c {criterion}')
        os.system('python compare_dfs.py')

if __name__ == '__main__':
    main()
