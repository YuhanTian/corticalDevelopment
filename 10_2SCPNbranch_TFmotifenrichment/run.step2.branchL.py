#!/usr/bin/env python
# coding=utf-8

# conda create -n scikit-learn
# conda activate scikit-learn
# conda install -c anaconda scikit-learn

import os
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import GridSearchCV
from pandas import read_csv
import pandas as pd


f = open('node2branchLMatrx_cell_gene_class.csv',encoding = 'utf-8')
data = read_csv(f)

X_train,y_train = data.iloc[:,2:], data.iloc[:,1]
y_train[~(y_train == "1early")] = "Other"

parameters = {'max_depth':[3,4,5], 'n_estimators':[25,50,75,100]}
clf = GradientBoostingClassifier(learning_rate = 1)
clf = GridSearchCV(clf, parameters, cv = 5, verbose = 5)
clf.fit(X_train,y_train)

print('best parameter of GBC:\n', clf.best_params_)

clf_best = GradientBoostingClassifier(n_estimators = 100, learning_rate = 1.0, max_depth = 5, random_state = 0).fit(X_train, y_train)
clf_best.score(X_train, y_train)
print('feature importances of GBC:\n',clf_best.feature_importances_)

feature_importance_table = pd.DataFrame(clf_best.feature_importances_, index = X_train.columns, columns=["importance"])
feature_importance_table.to_csv("node2_branchL_featuresImportance.txt", sep='\t')
