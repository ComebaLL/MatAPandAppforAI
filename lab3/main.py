import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.tsa.api import VAR
from statsmodels.tsa.stattools import adfuller


CSV_PATH = "energy_dataset.csv"
COLS = ["generation wind onshore", "total load actual", "price actual"]
FREQ = "D"          # частота агрегации: "D" - сутки 
MAX_LAGS = 15       # верхняя граница поиска p по информационному критерию
D_MAX = 1           # максимальный порядок интеграции; None -> оценить по ADF
ALPHA = 0.05        # уровень значимости

# OLS - метод наименьших квадратов
# VAR - векторная авторегрессия
# AIC - информативный критерий Акаике
# ADF - расширенный тест Дики-ФУллера

def load_data(path, cols, freq):
    """Читает CSV, оставляет нужные столбцы"""
    df = pd.read_csv(path, parse_dates=["time"])
    df["time"] = pd.to_datetime(df["time"], utc=True)
    df = df.set_index("time")[cols]
    return df.interpolate(limit_direction="both").resample(freq).mean().dropna()


def estimate_d_max(data, maxd=3, alpha=ALPHA):
    """Оценивает максимальный порядок интеграции по ADF.

    Ряд дифференцируется, пока ADF не отвергнет единичный корень. Сама модель
    дифференцированные ряды не использует: результат нужен только как число
    добавочных лагов. Это ориентир
    """
    def order(s):
        for d in range(maxd + 1):
            if adfuller(s.dropna(), autolag="AIC")[1] < alpha:
                return d            # p-value мал -> единичного корня нет
            s = s.diff()            # иначе берём первую разность и проверяем снова
        return maxd
    orders = {c: order(data[c]) for c in data.columns}
    return max(orders.values()), orders


def toda_yamamoto(data, caused, causing, p, d_max, robust=False):
    """Тест Тоды-Ямамото

    d_max - число добавочных лагов, которые в модели есть, но в гипотезу не входят
    causing - столбцы-следствия
    caused - столбцы-причины
    Алгоритм:
      1. Оцениваем уравнение для caused по OLS на уровнях, с k = p + d_max
         лагами всех переменных системы.
      2. Тестируем Вальдом, что коэффициенты при первых p лагах causing равны
         нулю. Добавочные d_max лагов остаются в модели, но в гипотезу не входят:
         именно они обеспечивают X^2(p)-распределение статистики при нестационарных
         рядах.
    H0: causing не является причиной по Грейнджеру для caused.
    """
    k = p + d_max

    # Матрица регрессоров: лаги 1..k каждой переменной + константа
    X = pd.concat({f"{c}.L{l}": data[c].shift(l)
                   for c in data.columns for l in range(1, k + 1)}, axis=1)
    X = sm.add_constant(X)

    # Склеиваем с зависимой переменной и убираем первые k строк с NaN от сдвигов
    m = pd.concat([data[caused], X], axis=1).dropna()

    model = sm.OLS(m[caused], m.drop(columns=caused))
    # HAC-ошибки (Ньюи-Уэста) защищают от гетероскедастичности и автокорреляции
    # остатков; для цен это актуально
    res = model.fit(cov_type="HAC", cov_kwds={"maxlags": k}) if robust else model.fit()

    # Матрица ограничений R*b = 0: по одной строке на каждый из первых p лагов
    names = list(res.params.index)
    R = np.zeros((p, len(names)))
    for i, l in enumerate(range(1, p + 1)):
        R[i, names.index(f"{causing}.L{l}")] = 1

    # use_f=False -> статистика X^2, число степеней свободы равно p
    w = res.wald_test(R, use_f=False, scalar=True)
    return float(w.statistic), float(w.pvalue)


def main():
    #  Данные 
    data = load_data(CSV_PATH, COLS, FREQ)
    print(f"Наблюдений: {len(data)}, период: {data.index[0].date()} – {data.index[-1].date()}")

    # Число добавочных лагов d_max 
    if D_MAX is None:
        d_max, orders = estimate_d_max(data)
        print("Порядки интеграции (ADF):", orders)
    else:
        d_max = D_MAX
    print(f"d_max = {d_max}")

    # Порядок VAR p по AIC на уровнях 
    p = VAR(data).select_order(maxlags=MAX_LAGS).selected_orders["aic"]
    p = max(p, d_max, 1)
    print(f"p (AIC) = {p}, лагов в VAR: {p + d_max}\n")

    # Тест для каждой упорядоченной пары переменных
    rows = []
    for caused in COLS:
        for causing in COLS:
            if caused == causing:
                continue
            chi2, pv = toda_yamamoto(data, caused, causing, p, d_max)
            rows.append({
                "гипотеза": f"{causing} -> {caused}",
                "chi2": round(chi2, 2),
                "df": p,                                   
                "p-value": round(pv, 4),
                "причинность": "да" if pv < ALPHA else "нет",
            })


    result = pd.DataFrame(rows)
    print(result.to_string(index=False))
    result.to_csv("toda_yamamoto_results.csv", index=False)


if __name__ == "__main__":
    main()