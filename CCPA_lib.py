

import itertools
import sys, os

from scipy.special import comb
from scipy import stats
import scipy.cluster.hierarchy as hac
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from matplotlib.colors import ListedColormap
from mpl_toolkits.mplot3d import Axes3D
from sklearn.linear_model import LinearRegression
from sklearn import metrics
from sklearn.ensemble import RandomForestClassifier
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches


def load_experiment_csvs(data_dpath=None, csv_fnames=None, meta_fnames=None):
    if data_dpath is None:
        data_dpath = r'C:\Users\wosnat\Documents\msc\work\CCPA\data_5_3_2019'
    if csv_fnames is None:
        # list of csvs - order is important
        csv_fnames = [
             'exp_1.csv',
             'exp2.csv',
             'exp3.csv',
             'exp4.csv',
             'exp5.csv',
             'exp6.csv',
             'exp7.csv',
        ]
    if meta_fnames is None:
        meta_fnames = [
             'exp_1_metadata.csv',
             'exp_2_metadata.csv',
             'exp_3_metadata.csv',
             'exp_4_metadata.csv',
             'exp_5_metadata.csv',
             'exp_6_metadata.csv',
             'exp_7_metadata.csv',
         ]

    dfs = [pd.read_csv(os.path.join(data_dpath, f)) for f in csv_fnames]
    meta_dfs = [pd.read_csv(os.path.join(data_dpath, f)) for f in meta_fnames]
    for i, d in enumerate(dfs, 1):
        d.rename(columns={'Unnamed: 0':'day'}, inplace=True)
        d['experiment'] = f'e{i}'

    melted_dfs = [d.melt(id_vars=['day', 'experiment'], var_name='sample', value_name='FL') for d in dfs]
    merged_dfs = [pd.merge(d, m, on='sample') for d,m in zip(melted_dfs, meta_dfs)]
    df = pd.concat(merged_dfs)
    return df

def update_calculated_fields(df, group_col=None, min_fl=0.05):
    if group_col is None:
        group_col = ['experiment', 'sample']
    df.loc[:, 'FL_orig'] = df.FL
    if min_fl is not None:
        df.loc[:, 'FL'] = df.FL.clip(lower=min_fl)

    df.loc[:, 'logFL'] =np.log(df['FL'])
    df.loc[:, 'cumsumFL'] = df.groupby(group_col).FL.transform(pd.Series.cumsum)
    df.loc[:, 'cumsumlogFL'] = df.groupby(group_col).logFL.transform(pd.Series.cumsum)
    df.loc[:, 'zscoreFL'] = df.groupby(group_col).FL.transform(stats.zscore)
    df.loc[:, 'diffFL'] = df.groupby(group_col).FL.transform(pd.Series.diff)
    df.loc[:, 'difflogFL'] = df.groupby(group_col).logFL.transform(pd.Series.diff)
    df.loc[:, 'diffday'] = df.groupby(group_col).day.transform(pd.Series.diff)
    df.loc[:, 'rateFL'] = df.diffFL / df.diffday
    df.loc[:, 'ratelogFL'] = df.difflogFL / df.diffday
    df.loc[:, 'experiment_sample'] = df.experiment + ', ' + df['sample']

    return df


def update_day_column(d, prev_exp, cur_exp):
    emax = d.loc[d.experiment == prev_exp].day.max()
    d.loc[d.experiment == cur_exp, 'day'] = d.loc[d.experiment == cur_exp, 'day'] + emax + 1
    return d

def concat_experiments(df):
    dfcat = df.loc[(df.experiment == 'e1') |
                   (df.experiment == 'e3') |
                   (df.experiment == 'e4') |
                   (df.experiment == 'e5') |
                   (df.experiment == 'e6')
                  ]
    dfcat = update_day_column(dfcat, 'e1', 'e3')
    dfcat = update_day_column(dfcat, 'e3', 'e4')
    dfcat = update_day_column(dfcat, 'e4', 'e5')
    dfcat = update_day_column(dfcat, 'e5', 'e6')
    #df.loc[:,'experiment'] = 'all'

    dfcat = update_calculated_fields(dfcat, group_col='sample')
    return dfcat

# Compute the following features based on a growth curve:
# 
# 1. lag phase: time until the first growth
#    - lag phase end time
#    - lag phase max value
#    (startThreashold = Application.WorksheetFunction.Max(yvalues) + 4 * (Application.WorksheetFunction.StDev(yvalues)))
# 2. growth phase: growth till max value
#    - growth max value
#    - growth diff (growth max value - lag phase max value)
#    - growth end time
#    - growth time (growth end time - lag phase end time
#    - growth rate (growth diff / growth time)
#    
# 3. stationary phase
#    -
# 
# 
# 4. decay phase
#    - decay time - total time to get to 0.5
#    - half yield - time to get to 1/2 growth max value
#        (halfdecayTime = (halfYield - beforeHalf + (decayPointsSlope * beforeHalfTime)) / decayPointsSlope)
#        (halfLife = halfdecayTime - maxTime)
#    - tenth yield - time to get to 1/10 growth max value
#    - decay rate
#    - number of spikes
#    
# 
# 5. global params
#   - AUC - area under the curve

def find_in_df(df, x_col, y_col, value):
    """ return time, value for the value closest to value in the value_col column
    """
    idx = df[y_col].sub(value).abs().idxmin()
    return df.loc[idx][[x_col, y_col]].tolist()

def find_midpoint(df, max_val, min_val, x_col, y_col, relative_coef):
    midpoint = (max_val - min_val)*relative_coef + min_val
    first_x = df[x_col].values[0]
    x, y = find_in_df(df, x_col, y_col, midpoint)
    return (x - first_x, y)
    
def fit_regression(df, x_col, y_col):
    """ fit linear regression. 
    return: model, intercept, coef, score
    """
    X = df[x_col].values.reshape((df.shape[0], 1))
    y = df[y_col]
    reg = LinearRegression(normalize=True).fit(X, y)
    return (reg, reg.intercept_, reg.coef_[0], reg.score(X, y))


def analyze_half_curve(df, label, max_val, min_val, x_col, y_col, y_logcol, rate_col, rate_logcol):
    half_day,     half_val    = find_midpoint(df, max_val, min_val, x_col, y_col, 0.5)
    loghalf_day,  loghalf_val = find_midpoint(df, np.log(max_val), np.log(min_val), x_col, y_logcol, 0.5)
    tenth_day,    tenth_val   = find_midpoint(df, max_val, min_val, x_col, y_col, 0.1)
    logtenth_day, logtenth_val = find_midpoint(df, np.log(max_val), np.log(min_val), x_col, y_logcol, 0.1)
    reg, intercept, coef, r2_score = fit_regression(df, x_col=x_col, y_col=y_logcol)
    try:
        auc = metrics.auc(df[x_col], df[y_col])
    except ValueError:
        auc = 0
    try:
        logauc = metrics.auc(df[x_col], df[y_logcol])
    except ValueError:
        logauc = 0
    diffval = df[y_col].max() - df[y_col].min()
    difflogval = df[y_logcol].max() - df[y_logcol].min()
    
    res = {
        f'{label}_half_day':     half_day,
        f'{label}_half_val':     half_val,
        f'{label}_log_half_day':  loghalf_day,
        f'{label}_log_half_val':  loghalf_val,
        f'{label}_tenth_day':    tenth_day,
        f'{label}_tenth_val':    tenth_val,
        f'{label}_log_tenth_day': logtenth_day,
        f'{label}_log_tenth_val': logtenth_val,
        f'{label}_log_intercept': intercept,
        f'{label}_log_coefficient': coef,
        f'{label}_log_score_r2': r2_score,
        f'{label}_auc': auc,
        f'{label}_logauc': logauc,
        f'{label}_diff': diffval,
        f'{label}_diff_log': difflogval,
        
        f'{label}_mean': df[y_col].mean(),
        f'{label}_std':  df[y_col].std(),
        f'{label}_min':  df[y_col].min(),
        f'{label}_max':  df[y_col].max(),
        
        f'{label}_mean_log': df[y_logcol].mean(),
        f'{label}_std_log':  df[y_logcol].std(),
        f'{label}_min_log':  df[y_logcol].min(),
        f'{label}_max_log':  df[y_logcol].max(),
        
        f'{label}_mean_rate': df[rate_col].mean(),
        f'{label}_std_rate':  df[rate_col].std(),
        f'{label}_min_rate':  df[rate_col].min(),
        f'{label}_max_rate':  df[rate_col].max(),
        
        f'{label}_mean_log_rate': df[rate_logcol].mean(),
        f'{label}_std_log_rate':  df[rate_logcol].std(),
        f'{label}_min_log_rate':  df[rate_logcol].min(),
        f'{label}_max_log_rate':  df[rate_logcol].max(),        
    }
    return res

def analyze_curve(df, x_col='day', y_col='FL', y_logcol='logFL', rate_col='rateFL', rate_logcol='ratelogFL', meta_col=None):
    if meta_col is None:
        meta_col = ['experiment', 'sample', 'PRO', 'ALT', 'culture']
    df = df.reset_index()
    maxval = df[y_col].max()
    minval = df[y_col].min()
    maxlogval = df[y_logcol].max()
    minlogval = df[y_logcol].min()
    maxidx = df[y_col].idxmax()
    maxday = df.loc[maxidx][x_col]
    res = {
        'max' : maxval,
        'max_log' : maxlogval,
        'max_day' : maxday,
        'min' : minval,
        'min_log' : minlogval,
    }
    for c in meta_col:
        res[c] = df[c].unique()[0]
    res.update(analyze_half_curve(df[:maxidx+1], label='growth', 
                                  max_val=maxval, min_val=minval,
                                  x_col=x_col, y_col=y_col, y_logcol=y_logcol,
                                  rate_col=rate_col, rate_logcol=rate_logcol,
                                 ))
    res.update(analyze_half_curve(df[maxidx:], label='decay',
                                  max_val=maxval, min_val=minval,
                                  x_col=x_col, y_col=y_col, y_logcol=y_logcol,
                                  rate_col=rate_col, rate_logcol=rate_logcol,
                                 ))
    return pd.Series(res)


def display_features(d, res, x_col='day', y_col='FL', y_logcol='logFL'):
    sns.lineplot(data=d, x=x_col, y=y_col, label='FL')
    sns.lineplot(data=d, x=x_col, y=y_logcol, label='logFL')
    sns.scatterplot(x=[res['max_day']], y=[res['max']], label='max', size=[100], legend=False)
    sns.scatterplot(x=[res['max_day']], y=[res['max_log']], label='maxlog', size=[100], legend=False)

    sns.scatterplot(x=[res['growth_half_day']], y=[res['growth_half_val']],  size=[100], legend=False)
    sns.scatterplot(x=[res['growth_log_half_day']], y=[res['growth_log_half_val']],  size=[100], legend=False)
    sns.scatterplot(x=[res['growth_tenth_day']], y=[res['growth_tenth_val']],  size=[100], legend=False)
    sns.scatterplot(x=[res['growth_log_tenth_day']], y=[res['growth_log_tenth_val']],  size=[100], legend=False)
    up_x = d.loc[d[x_col] <= res['max_day'], x_col]
    up_pred = up_x * res['growth_log_coefficient'] + res['growth_log_intercept']
    sns.lineplot(x=up_x, y=up_pred)

    sns.scatterplot(x=[res['max_day']+res['decay_half_day']], y=[res['decay_half_val']],  size=[100], legend=False)
    sns.scatterplot(x=[res['max_day']+res['decay_log_half_day']], y=[res['decay_log_half_val']],  size=[100], legend=False)
    sns.scatterplot(x=[res['max_day']+res['decay_tenth_day']], y=[res['decay_tenth_val']],  size=[100], legend=False)
    sns.scatterplot(x=[res['max_day']+res['decay_log_tenth_day']], y=[res['decay_log_tenth_val']], size=[100], legend=False)
    dn_x = d.loc[d[x_col] >= res['max_day'], x_col]
    dn_pred = dn_x * res['decay_log_coefficient'] + res['decay_log_intercept']
    sns.lineplot(x=dn_x, y=dn_pred)


def generate_features(df, groupby_col='experiment_sample'):
    g = df.groupby(groupby_col).apply(analyze_curve)
    return g

def run_pca(X, metadf, sample_col='experiment_sample', n_components=2):
    scaledX = StandardScaler().fit_transform(X)
    pca = PCA(n_components=n_components)
    principalComponents = pca.fit_transform(scaledX)
    print('Variance percent explained\n', pca.explained_variance_ratio_)
    pca_columns = [f'PCA{i}' for i in range(1,n_components+1)]
    principalDf = pd.DataFrame(data = principalComponents
                 , columns = pca_columns)
    principalDf.set_index(X.index, inplace=True)
    dfpca = pd.merge(left=principalDf, left_index=True, right=metadf, right_on=sample_col)
    return dfpca


def display_pca_3d(dfpca, color_col, style_col):
    fig = plt.figure(figsize=(15, 10))
    ax = fig.add_subplot(111, projection='3d')
    numcolors = dfpca[color_col].nunique()
    color_labels = dfpca[color_col].unique().tolist()
    numshapes = dfpca[style_col].nunique()
    style_labels = dfpca[style_col].unique().tolist()
    markers = Line2D.filled_markers
    pallete = sns.color_palette('hls', n_colors=numcolors).as_hex()
    for i, a in enumerate(color_labels):
        for j, b in enumerate(style_labels):
            d = dfpca.loc[(dfpca[color_col] == a) & (dfpca[style_col] == b)]
            ax.scatter(d.PCA1, d.PCA2, d.PCA3, s=300, alpha=0.4, edgecolors='w', c=pallete[i], marker=markers[j],)
    rows = [mpatches.Patch(color=pallete[i]) for i in range(numcolors)]
    columns = [plt.plot([], [], markers[i], markerfacecolor='w',
                        markeredgecolor='k')[0] for i in range(numshapes)]
    plt.legend(rows + columns, color_labels + style_labels, loc=2, bbox_to_anchor=(1.05, 1), borderaxespad=0.)
    #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.set_zticklabels([])
    return ax

def get_meta(df, meta_col=None, value_col='FL'):
    if meta_col is None:
        meta_col = ['experiment_sample', 'experiment', 'sample', 'PRO', 'ALT', 'culture']
    metadf = df.groupby(meta_col)[value_col].mean().reset_index()
    # todo meta.drop value_col
    return metadf


def dfcat2X(df, sample_col='sample', value_col='cumsumlogFL', x_col='day' ):
    X = df.pivot(index=sample_col, columns=x_col, values=value_col)
    X.fillna(method='pad', inplace=True)
    return X

    # # todo meta
    # # meta = meta_dfs[0]
    # # dfpca = pd.merge(left=principalDf, left_index=True, right=meta, right_on='sample')

def features2X(df, meta_cols=None):
    if meta_cols is None:
        meta_cols = ['sample', 'experiment','PRO', 'ALT', 'culture', 'experiment_sample']

    feature_cols = df.columns.difference(meta_cols)
    X = df[feature_cols]
    X.fillna(value=0, inplace=True)
    return X

def experiments2X(df, sample_col='experiment_sample', value_col='logFL', x_col='day', cumsummode=True):
    X = df.pivot(index=sample_col, columns=x_col, values=value_col)
    Xi = X.interpolate(method='from_derivatives', axis=1, limit_area='inside')
    Xi.fillna(method='pad', inplace=True, axis=1)
    if cumsummode:
        Xc = Xi.transform(pd.Series.cumsum, axis=1)
    else:
        Xc = Xi
    return Xc


def add_metacols_to_pca(dfpca, df, meta_cols):
    for c in meta_cols:
        dfpca[c] = df[c]
    return dfpca


def forest_classifier(X, y):
    scaledX = StandardScaler().fit_transform(X)
    clf = RandomForestClassifier(n_estimators=100, max_depth=2,
                                  random_state=0)
    clf.fit(scaledX, y)
    print(clf.score(scaledX, y))
    return clf

def forest_feature_importance(clf, col_names, n=10):
    feature_importances = pd.DataFrame(clf.feature_importances_,
                                       index = col_names,
                                        columns=['importance']).sort_values('importance', ascending=False)
    feature_importances.nlargest(columns='importance',n=n).plot(kind='barh')


def forest_heatmap(clf, X, y, metadf=None, breakdown=None):
    scaledX = StandardScaler().fit_transform(X)
    d = pd.DataFrame()
    d['actual'] = y
    d['predicted'] = clf.predict(scaledX)
    d['x'] =1
    idx = 'actual'
    if metadf is not None and breakdown is not None:
        d[breakdown] = metadf[breakdown]
        idx = [breakdown, 'actual']

    t =d.pivot_table(index=idx,columns=['predicted'], aggfunc='count')
    t.columns = t.columns.get_level_values(1)
    sns.heatmap(t,annot=True)


def extract_decay(d, y_col = 'FL', x_col = 'day', scale=True):
    y_col_orig = f'{y_col}_orig'
    d = d.reset_index()
    maxidx = d[y_col].idxmax()
    c = d[maxidx:].reset_index()
    r = pd.DataFrame()

    r[x_col] = c[x_col] - c[x_col].min()
    if scale:
        r[y_col] = c[y_col] / c[y_col].max()
    else:
        r[y_col] = c[y_col]
    
    meta_cols = ['experiment', 'sample', 'PRO', 'ALT', 'culture']
    for c in meta_cols:
        r[c] = d[c].unique()[0]
    return r


def generate_decay(df, sample_col='experiment_sample', scale=True):
    dfd = df.groupby(sample_col).apply(lambda x: extract_decay(x, scale=scale))
    dfd.reset_index(level=0, inplace=True)
    return dfd



if __name__ == '__main__':

    df = pd.read_pickle('CCPA.pkl.gz')


    d = df.loc[(df.experiment == 'e5') & (df['sample'] == '37C')]
    analyze_curve(d)

