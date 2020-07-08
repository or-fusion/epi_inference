import sys
import pandas as pd
import plotnine as p9


def ribbon_plot(fname, title, infectious_period=None):
    df = pd.read_csv(fname, parse_dates=['dates'])
    fvi = df['qmean_filtered_est_beta'].first_valid_index()
    if fvi is None:
        print('No data found in {}'.format(fname))
        return

    df = df.loc[fvi:]
    df['plot_mean'] = df['qmean_filtered_est_beta'] * infectious_period
    df['plot_low'] = df['q05_filtered_est_beta'] * infectious_period
    df['plot_high'] = df['q95_filtered_est_beta'] * infectious_period

    myplot = (p9.ggplot(data=df, mapping=p9.aes(x='dates', y='plot_mean')) +
              p9.geom_line() +
              p9.geom_ribbon(data=df, mapping=p9.aes(ymin='plot_low',
                                                     ymax='plot_high'), alpha=0.3) +
              p9.coord_cartesian(ylim=(0, 4)) +
              p9.xlab('') +
              p9.ylab('Estimated Transmission Rate') +
              p9.ggtitle('County FIPS: ' + str(title)))
    return myplot


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: python plot_beta_summary.py <summary_estimated_beta_county_?????.csv> <plot title>')
        quit()

    fname = sys.argv[1]
    print(fname)
    title = sys.argv[1][-9:-4]
    plt = ribbon_plot(fname, title, infectious_period=4.3)
    print(plt)