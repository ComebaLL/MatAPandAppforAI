import numpy as np
import matplotlib
import matplotlib.pyplot as plt

UNVISITED = -2   # точка ещё не обработана
NOISE = -1       # шум


def region_query(X, i, eps):
    """eps-окрестность точки i: индексы всех q, для которых d(X[i], q) <= eps.
    Сама точка входит в свою окрестность"""
    dist = np.sqrt(np.sum((X - X[i]) ** 2, axis=1))  # евклидово расстояние
    return np.where(dist <= eps)[0]


def dbscan(X, eps, min_samples):
    """Простейший DBSCAN.

    X           -- массив формы (n, d)
    eps         -- радиус окрестности
    min_samples -- минимум точек в eps-окрестности для ядровой точки
    Возвращает массив меток: 0, 1, 2, ... -- кластеры, -1 -- шум.
    """
    X = np.asarray(X, dtype=float)
    n = len(X)
    labels = np.full(n, UNVISITED)
    cluster_id = 0

    for p in range(n):
        if labels[p] != UNVISITED:          # уже обработана -- пропускаем
            continue

        neighbors = region_query(X, p, eps)

        if len(neighbors) < min_samples:    # не ядровая -> пока шум
            labels[p] = NOISE               # (может стать граничной позже)
            continue

        # p -- ядровая точка: создаём новый кластер и расширяем его (BFS)
        labels[p] = cluster_id
        queue = list(neighbors)
        k = 0
        while k < len(queue):
            q = queue[k]
            k += 1

            if labels[q] == NOISE:          # шум, достижимый из ядровой,
                labels[q] = cluster_id      # -> граничная точка кластера
                continue
            if labels[q] != UNVISITED:      # уже принадлежит кластеру
                continue

            labels[q] = cluster_id
            q_neighbors = region_query(X, q, eps)
            if len(q_neighbors) >= min_samples:   # q ядровая -> расширяем
                queue.extend(q_neighbors)

        cluster_id += 1

    return labels


def k_distances(X, k):
    """Для подбора eps: отсортированные расстояния до k-го ближайшего соседа"""
    X = np.asarray(X, dtype=float)
    result = np.empty(len(X))
    for i in range(len(X)):
        dist = np.sqrt(np.sum((X - X[i]) ** 2, axis=1))
        result[i] = np.sort(dist)[k]  # [0] -- сама точка, поэтому индекс k
    return np.sort(result)

def plot_k_distance(X, k, ax):
    """График k-расстояний: «локоть» подсказывает значение eps."""
    kd = k_distances(X, k)
    ax.plot(kd)
    ax.set_title(f"k-distance graph (k = {k})")
    ax.set_xlabel("точки, отсортированные по k-расстоянию")
    ax.set_ylabel(f"расстояние до {k}-го соседа")
    ax.grid(alpha=0.3)
 
 
def plot_clusters(X, labels, ax, eps, min_samples):
    """Точки раскрашены по кластерам, шум -- серые крестики."""
    for c in sorted(set(labels)):
        mask = labels == c
        if c == NOISE:
            ax.scatter(X[mask, 0], X[mask, 1], c="gray", marker="x",
                       label="шум")
        else:
            ax.scatter(X[mask, 0], X[mask, 1], s=30, label=f"кластер {c}")
    ax.set_title(f"DBSCAN: eps = {eps}, min_samples = {min_samples}")
    ax.set_xlabel("x1")
    ax.set_ylabel("x2")
    ax.legend()
    ax.grid(alpha=0.3)


if __name__ == "__main__":
    rng = np.random.default_rng(42)

    # два облака + немного равномерного шума
    cluster1 = rng.normal(loc=[0, 0], scale=0.3, size=(50, 2))
    cluster2 = rng.normal(loc=[3, 3], scale=0.3, size=(50, 2))
    noise = rng.uniform(low=-2, high=5, size=(10, 2))
    X = np.vstack([cluster1, cluster2, noise])

    labels = dbscan(X, eps=0.4, min_samples=4)

    print("Кластеры:", sorted(set(labels) - {NOISE}))
    print("Размеры кластеров:",
          {int(c): int(np.sum(labels == c)) for c in set(labels) if c != NOISE})
    print("Точек шума:", int(np.sum(labels == NOISE)))
    print("k-distance (каждая 20-я):", np.round(k_distances(X, 4)[::20], 3))

    # графики слева k-distance, справа результат кластеризации
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    plot_k_distance(X, 4, axes[0])
    axes[0].axhline(0.4, color="red", linestyle="--", label="выбранный eps")
    axes[0].legend()
    plot_clusters(X, labels, axes[1], eps=0.4, min_samples=4)
    plt.tight_layout()
    plt.savefig("dbscan_result.png", dpi=120)
    plt.show()