sudo sed  -i '' 's/pd.tslib/pd/g'  /Library/Python/3.7/site-packages/ggplot/utils.py
sudo sed  -i '' 's/pd.tslib/pd/g'  /Library/Python/3.7/site-packages/ggplot/stats/smoothers.py
sudo sed  -i '' 's/from pandas.lib import Timestamp/from pandas import Timestamp/g'  /Library/Python/3.7/site-packages/ggplot/stats/smoothers.py

