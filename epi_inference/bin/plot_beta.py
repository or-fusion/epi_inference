import sys
import pandas as pd
import matplotlib.pyplot as plt

def add_line_ribbon(ax, fname, ribbon=False, infectious_period=None):
    df = pd.read_csv(fname, parse_dates=['dates'])
    df.dropna(axis=0, inplace=True)    
    x = df['dates'].values
    ym = df['qmean_filtered_est_beta'].values*infectious_period
    if ribbon:
        y05 = df['q05_filtered_est_beta'].values*infectious_period
        y95 = df['q95_filtered_est_beta']*infectious_period
        ax.fill_between(x, y05, y95, interpolate=True, color='silver')

    ax.plot(x,ym)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: python plot_beta_summary.py <summary_estimated_beta_county_?????.csv> <plot title>')
        quit()

    fig, ax = plt.subplots(1, 1, sharex=True)
    plt.subplots_adjust(bottom=0.2)

    fname = sys.argv[1]
    print(fname)
    add_line_ribbon(ax, fname, ribbon=True, infectious_period=4.3)
    for fname in sys.argv[2:-1]:
        print(fname)
        add_line_ribbon(ax, fname, ribbon=False, infectious_period=4.3)

    title = sys.argv[-1]
    print(title)
    plt.title(title)
    plt.xticks(rotation=45)
    plt.show()



