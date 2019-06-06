

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
from sklearn.metrics import classification_report, accuracy_score

from scipy.optimize import least_squares
from scipy.optimize import curve_fit

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

def update_rolling_average(df, group_col=None, min_fl=0.05):
    if group_col is None:
        group_col = ['experiment', 'sample']
    df = df.groupby(group_col).apply(add_rolling_average)
    return df.reset_index(drop=True)


def add_rolling_average(df, from_col='FL', to_col='roll', index_col='day', window='3d'):
    """ add a rolling average column to df """
    #print(df.head())
    c = df.copy()
    #print(c.head())
    c.index = pd.to_timedelta(c[index_col], unit='d')
    c.index.names = ['i'] # change from 'day' to avoid merge conflict
    c.loc[:, to_col] = c[from_col].rolling(window, min_periods=1).mean()
    d = df.merge(c.loc[:, [index_col, to_col]], on=index_col, suffixes=('',''))
    #print(d[to_col])
    return d


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
    d = df.dropna(subset=[value_col])
    X = d.pivot(index=sample_col, columns=x_col, values=value_col)
    if X.columns.dtype == '<m8[ns]':
        X.columns = X.columns.astype('timedelta64[D]')
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


def _resample_func(df, x_col='day', value_col='FL', period='3d' ):
    t = df
    t.index = pd.to_timedelta(t[x_col], unit='d')
    return t.resample(period).agg({value_col : 'mean'})
    #return t.resample(period).agg({y_col : ['mean', 'median','std']})
    #return t.rolling(period, min_periods=1).agg({y_col : ['mean', 'median','std']})


def resample_df(df, x_col='day', value_col='FL', period='3d', groupby_cols=None):
    if groupby_cols is None:
        groupby_cols = ['experiment_sample', 'experiment', 'sample', 'PRO', 'ALT', 'culture']
    df_resampled = df.groupby(groupby_cols).apply(lambda x: _resample_func(x, value_col=value_col, x_col=x_col, period=period))
    df_resampled = df_resampled.reset_index()
    df_resampled.dropna(inplace=True)
    return df_resampled



def forest_classifier(X, y):
    scaledX = StandardScaler().fit_transform(X)
    clf = RandomForestClassifier(n_estimators=100, oob_score=True,
                                 #max_depth=2,
                                 # random_state=0)
                                 )
    clf.fit(scaledX, y)
    print('train score', clf.score(scaledX, y))
    print ('oob score', clf.oob_score_)
    return clf

def forest_feature_importance(clf, col_names, n=10):
    feature_importances = pd.DataFrame(clf.feature_importances_,
                                       index = col_names,
                                        columns=['importance']).sort_values('importance', ascending=False)
    feature_importances.nlargest(columns='importance',n=n).plot(kind='barh')


def forest_heatmap(clf, X, y, metadf=None, breakdown=None, func=None, ax=None):
    scaledX = StandardScaler().fit_transform(X)
    d = pd.DataFrame()
    d['actual'] = y
    d['predicted'] = clf.predict(scaledX)
    if func is not None:
        d['actual'] = func(d['actual'])
        d['predicted'] = func(d['predicted'])

    d['x'] = 1
    idx = 'actual'
    if metadf is not None and breakdown is not None:
        d[breakdown] = metadf[breakdown]
        idx = [breakdown, 'actual']

    t = d.pivot_table(index=idx, columns=['predicted'], aggfunc='count')
    t.columns = t.columns.get_level_values(1)
    print(f"accuracy: {accuracy_score(y_true=d['actual'], y_pred=d['predicted'])}")
    print(classification_report(y_true=d['actual'], y_pred=d['predicted']))

    #forest_feature_importance(clf, X.columns)

    sns.heatmap(t, annot=True, cmap='coolwarm', ax=ax)
    #return t

# #func = lambda x : x.str.split(', ', expand=True)[0]
# t =forest_heatmap(clf=clf, X=X, y=y)
# fig, ax = plt.subplots(figsize=(12,12))
# sns.heatmap(t, annot=True, cmap='coolwarm', ax=ax)
#
# ax.axhline(5)
# ax.axhline(15)
# ax.axhline(20)
# ax.axhline(10)
# ax.axvline(4)
# ax.axvline(8)
# ax.axvline(13)
# ax.axvline(15)


import sklearn.metrics as metrics


def _calc_score_for_one_type(res, y_true, y_pred, suffix):
    precision, recall, f1, support = metrics.precision_recall_fscore_support(
        y_true=y_true,
        y_pred=y_pred,
        average='weighted'
    )
    res[f'accuracy_{suffix}'] = metrics.accuracy_score(y_true=y_true, y_pred=y_pred)
    res[f'precision_{suffix}'] = precision
    res[f'recall_{suffix}'] = recall
    res[f'f1_{suffix}'] = f1
    res[f'support_{suffix}'] = support


def score_model(modelname, clf, X_train, y_train, X_test, y_test, return_y = False):
    func = lambda x: x.str.split(',', expand=True)[0]
    scalar = StandardScaler()
    scaledX_train = scalar.fit_transform(X_train)
    scaledX_test = scalar.transform(X_test)

    y_train_pred = clf.predict(scaledX_train)
    y_test_pred = clf.predict(scaledX_test)

    y_train_pro = func(y_train)
    y_test_pro = func(y_test)
    y_train_pred_pro = func(pd.Series(y_train_pred, index=y_train.index))
    y_test_pred_pro = func(pd.Series(y_test_pred, index=y_test.index))

    res = {'model': modelname}
    res['oob_score'] = clf.oob_score_
    _calc_score_for_one_type(res, y_true=y_train, y_pred=y_train_pred, suffix='train')
    _calc_score_for_one_type(res, y_true=y_test, y_pred=y_test_pred, suffix='test')
    _calc_score_for_one_type(res, y_true=y_train_pro, y_pred=y_train_pred_pro, suffix='train_PRO')
    _calc_score_for_one_type(res, y_true=y_test_pro, y_pred=y_test_pred_pro, suffix='test_PRO')

    if not return_y:
        return clf, res, None
    else:
        y_train_df = pd.DataFrame(data={
            f'{modelname}_y': y_train,
            f'{modelname}_y_PRO': y_train_pro,
            f'{modelname}_y_pred': y_train_pred,
            f'{modelname}_y_pred_PRO': y_train_pred_pro,
        }, index=y_train.index)
        y_train_df['Type'] = 'Train'
        y_test_df = pd.DataFrame(data={
            f'{modelname}_y': y_test,
            f'{modelname}_y_PRO': y_test_pro,
            f'{modelname}_y_pred': y_test_pred,
            f'{modelname}_y_pred_PRO': y_test_pred_pro,
        }, index=y_test.index)
        y_test_df['Type'] = 'Test'
        y_df = pd.concat([y_train_df, y_test_df])
        return clf, res, y_df


def ml(X, metadf, modelname, y_col='PRO'):
    print (modelname)
    X_train = X[X.index.str.startswith('e1') |
                X.index.str.startswith('e3') |
                X.index.str.startswith('e4') |
                X.index.str.startswith('e5')]
    X_test = X[X.index.str.startswith('e6')]

    metadf.index = metadf.experiment_sample
    metadf_train = metadf[metadf.index.str.startswith('e1') |
                          metadf.index.str.startswith('e3') |
                          metadf.index.str.startswith('e4') |
                          metadf.index.str.startswith('e5')]
    metadf_test = metadf[metadf.index.str.startswith('e6')]

    if y_col=='PRO_ALT':
        y_train = metadf_train.PRO + ',' + metadf_train.ALT
        y_test = metadf_test.PRO + ',' + metadf_test.ALT
    else:
        y_train = metadf_train[y_col]
        y_test = metadf_test[y_col]

    clf = forest_classifier(X=X_train, y=y_train)
    return score_model(modelname, clf, X_train, y_train, X_test, y_test)


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



#############################

# models
#
# J.J. Arps (1945):
#
# Exponential:  The decline rate a does not varies with q.
#
#               q = q0* e^(-a*t)
#
# Harmonic:     The decline rate a varies linearly with q.
#
#               q = q0 / (1 + a * t)
#
# Hyperbolic:   The decline rate a varies geometrically with q
#
#               q = q0 / (1 + d*a*t)^(1/d)
#
#
#  3-­‐parameter  Logistic :
#         Y(t)  =  b1   /  [1  +  b2 *exp(  b3 *t  )]
#
# 4-­‐parameter  Logistic
#
#         Y(t)  =  b1  +  (b2   –  b1 )  /  [1  +  exp(b3 *(t  –  b4 ))]
#
# 3-­‐parameter  Rodbard
#
# Y(t)  =  (1   –  b4 )/[1  +  (t/b3 ) b2 ]  +  b4
#
# 4-­‐parameter  Rodbard
#
#     Y(t)  =  (b1   –  b4 )/[1  +  (t/b3 ) b2 ]  +  b4
#
# Gompertz
#
#     Y(t)  =  b1 *exp(  –  b2* exp(  –  b3 *t))
#
# Log-­‐Logistic
#
#     Y(t)  =  b1  –  log(1  +  b2 *exp( –  b3 *t)
#
# First-­‐order  Decay
#
#     Y(t)  =  b1 *exp(  –  b2 *t)  +  b3
#
# s-curve
# g(x) = 1 - 1/(1 + exp(-x))

def model_exponential(z, a1, b1, c1, _):
    return (
        (a1 * np.exp(-b1 * z) + c1)
    )

def model_harmonic(z, a1, b1, c1, _):
    return (
        ((a1 / (1 + b1 * z)) + c1)
    )

def model_hyperbolic(z, a1, b1, c1, d1):
    return (
        ((a1 / np.power((1 + d1 * b1 * z), d1)) + c1)
    )

def model5(z,  # s1, #s2, s3,
           a1, b1, c1, d1  # a2, b2,c2, a3, b3, c3,  a4, b4, c4
           ):
    return (
        ((a1 * np.power((1 + d1 * b1 * z), -1 / d1)) + c1)
    )

def model_exponential_segmented(z, s1, s2, s3, a1, b1, c1, a2, b2, c2, a3, b3, c3, a4, b4, c4):
    return np.piecewise(z.values,
                        condlist=[
                            (z.values < s1),
                            (z.values >= s1) & (z.values < s2),
                            (z.values >= s2) & (z.values < s3),
                        ],
                        funclist=[
                            lambda t: model_exponential(t, a1, b1, c1,0),
                            lambda t: model_exponential(t, a2, b2, c2,0),
                            lambda t: model_exponential(t, a3, b3, c3,0),
                            lambda t: model_exponential(t, a4, b4, c4,0),
                        ]
                        )

def model_linear(z, a1, b1, _1, _2 ):
    return a1 * z + b1

#  3-­‐parameter  Logistic :
#         Y(t)  =  b1   /  [1  +  b2 *exp(  b3 *t  )]
def model_logistic3(z, b1, b2, b3 , _):
    return b1 / (1 + b2 * np.exp(b3 * z))
#
# 4-­‐parameter  Logistic
#
#         Y(t)  =  b1  +  (b2   –  b1 )  /  [1  +  exp(b3 *(t  –  b4 ))]
def model_logistic4(z, b1, b2, b3, b4 ):
    return b1 + (b2 - b1) / (1 +  np.exp(b3 * (z - b4)))
#
# 3-­‐parameter  Rodbard
#
# Y(t)  =  (1   –  b4 )/[1  +  (t/b3 ) b2 ]  +  b4
def model_rodbard4(z, b1, b2, b3, b4):
    return (b1 - b4) / (1 + (z / b3)* b2) + b4
#
# 4-­‐parameter  Rodbard
#
#     Y(t)  =  (b1   –  b4 )/[1  +  (t/b3 ) b2 ]  +  b4
#
# Gompertz
#
#     Y(t)  =  b1 *exp(  –  b2* exp(  –  b3 *t))
def model_gompertz(z, b1, b2, b3, _):
    return  b1 * np.exp( - b2 * np.exp( -b3 * z))
#
# Log-­‐Logistic
#
#     Y(t)  =  b1  –  log(1  +  b2 *exp( –  b3 *t)
def model_loglogistic(z, b1, b2, b3, _):
    return  b1 - np.log(1 + b2 * np.exp( -b3 * z))
#
# First-­‐order  Decay
#
#     Y(t)  =  b1 *exp(  –  b2 *t)  +  b3


def model_cubic(z, b1, b2, b3, b4):
    return  b1* np.power(z, 3) + b2 * np.power(z, 2) + b3 * z + b4

def model_scurve(z, b1, b2, b3, _):
    return b1 * (1 - 1 / (1 + np.exp(-b2 * z))) + b3


# if __name__ == '__main__':
#
#     df = pd.read_pickle('CCPA.pkl.gz')
#
#
#     d = df.loc[(df.experiment == 'e5') & (df['sample'] == '37C')]
#     analyze_curve(d)
#


# value_col_list = ['FL', 'logFL', 'cumsumFL', 'cumsumlogFL', 'zscoreFL', 'rateFL', 'ratelogFL'],
# resample_period_list = [None, '1d', '3d', '5d']
# y_col_list = ['PRO', 'ALT', 'PRO_ALT']
# cumsummode_list = [False, True]


def compare_models(df, value_col, resample_period_list, y_col_list, cumsummode_list):
    # resample first
    stats_list = []
    for resample_period in resample_period_list:
        if resample_period is not None:
            d = cp.resample_df(df, value_col=value_col, period=resample_period)
        else:
            d = df
        metadf = cp.get_meta(d, value_col=value_col)
        for cumsummode in cumsummode_list:
            X = cp.experiments2X(d, value_col=value_col, cumsummode=cumsummode)

            resample_str = '_' + resample_period if resample_period is not None else ''
            cumsummode_str = '_cumsum' if cumsummode else ''
            for y_col in y_col_list:
                modelname = f'{y_col}_{value_col}{resample_str}{cumsummode_str}'
                clf, res, _ = cp.ml(X=X, metadf=metadf, modelname=modelname, y_col=y_col)
                stats_list.append(res)
    return stats_list