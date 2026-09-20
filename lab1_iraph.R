# install.packages("igraph")

library(igraph)



# ============================================================
# 1. Создание случайного неориентированного графа
# ============================================================

set.seed(42)

n_vertices <- 20

# Полный граф K_n имеет n(n-1)/2 рёбер.
max_edges <- n_vertices * (n_vertices - 1) / 2

# Требуемая плотность примерно 1/3–1/2.
n_edges <- 70

graph_density_target <- n_edges / max_edges

# Генерируем граф до тех пор, пока он не станет связным.
repeat {
  g <- sample_gnm(
    n = n_vertices,
    m = n_edges,
    directed = FALSE,
    loops = FALSE
  )

  if (is_connected(g)) {
    break
  }
}

V(g)$name <- paste0("v", seq_len(vcount(g)))

cat("Количество вершин:", vcount(g), "\n")
cat("Количество рёбер:", ecount(g), "\n")
cat("Максимально возможное число рёбер:", max_edges, "\n")
cat("Плотность графа:", edge_density(g), "\n")
cat("Плотность относительно полного графа:", graph_density_target, "\n")
cat("Граф ориентированный:", is_directed(g), "\n")
cat("Граф связный:", is_connected(g), "\n")


# ============================================================
# 2. Визуализация графа
# ============================================================

set.seed(42)

plot(
  g,
  layout = layout_with_fr(g),
  vertex.label = V(g)$name,
  vertex.size = 18,
  vertex.color = "lightblue",
  vertex.frame.color = "black",
  edge.color = "grey50",
  main = "Случайный неориентированный граф"
)


# ============================================================
# 3. Связность графа
# ============================================================

cat("\n========== 1. ПРОВЕРКА СВЯЗНОСТИ ==========\n")

connected_result <- is_connected(g)

cat("Граф связный:", connected_result, "\n")

components_result <- components(g)

cat("Количество компонент связности:", components_result$no, "\n")
cat("Принадлежность вершин компонентам:\n")
print(components_result$membership)


# ============================================================
# 4. Проверка, является ли граф деревом
# ============================================================

cat("\n========== 2. ПРОВЕРКА НА ДЕРЕВО ==========\n")

tree_result <- is_tree(g)

cat("Граф является деревом:", tree_result, "\n")

# ============================================================
# 6. Кратчайшие пути
# ============================================================

cat("\n========== 4. КРАТЧАЙШИЕ ПУТИ ==========\n")

# Матрица длин кратчайших путей между всеми парами вершин.
distance_matrix <- distances(
  g,
  mode = "all",
  weights = NA
)

cat("Матрица расстояний между всеми вершинами:\n")
print(distance_matrix)

# Один конкретный кратчайший путь:
source_vertex <- 1
target_vertex <- n_vertices

shortest_path <- shortest_paths(
  g,
  from = source_vertex,
  to = target_vertex,
  mode = "all",
  weights = NA,
  output = "vpath"
)

cat(
  "\nОдин кратчайший путь от вершины",
  source_vertex,
  "до вершины",
  target_vertex,
  ":\n"
)

print(shortest_path$vpath[[1]])


# ============================================================
# 7. Соседи произвольной вершины
# ============================================================

cat("\n========== 5. СОСЕДИ ВЕРШИНЫ ==========\n")

selected_vertex <- 1

neighbors_result <- neighbors(
  g,
  v = selected_vertex,
  mode = "all"
)

cat(
  "Соседи вершины",
  selected_vertex,
  ":\n"
)

print(neighbors_result)


# ============================================================
# 8. Все треугольники, содержащие выбранную вершину
# ============================================================

cat("\n========== 6. ТРЕУГОЛЬНИКИ ==========\n")

selected_vertex <- 1

all_triangles <- triangles(g)

# triangles() возвращает вершины треугольников
# в виде последовательности по 3 вершины.
if (length(all_triangles) > 0) {

  triangle_matrix <- matrix(
    all_triangles,
    nrow = 3
  )

  selected_triangles <- triangle_matrix[
    ,
    apply(
      triangle_matrix,
      2,
      function(x) selected_vertex %in% x
    ),
    drop = FALSE
  ]

  cat(
    "Треугольники, содержащие вершину",
    selected_vertex,
    ":\n"
  )

  if (ncol(selected_triangles) > 0) {
    print(selected_triangles)
  } else {
    cat("У выбранной вершины нет треугольников.\n")
  }

} else {

  cat("В графе нет треугольников.\n")
}

cat(
  "\nКоличество треугольников для каждой вершины:\n"
)

print(count_triangles(g))


# ============================================================
# 9. Центральности вершин и рёбер
# ============================================================

cat("\n========== 7. ЦЕНТРАЛЬНОСТИ ==========\n")

centrality_table <- data.frame(
  vertex = seq_len(vcount(g)),
  degree = degree(g, mode = "all"),
  betweenness = betweenness(
    g,
    directed = FALSE,
    normalized = FALSE
  ),
  closeness = closeness(
    g,
    mode = "all",
    normalized = TRUE,
    weights = NA
  ),
  eigenvector = eigen_centrality(
    g,
    directed = FALSE,
    weights = NA
  )$vector,
  harmonic = harmonic_centrality(
    g,
    mode = "all",
    weights = NA
  ),
  page_rank = page_rank(
    g,
    directed = FALSE,
    weights = NA
  )$vector
)

# Дополнительные меры центральности.
alpha_result <- tryCatch(
  alpha_centrality(
    g,
    alpha = 0.01,
    weights = NA
  ),
  error = function(e) rep(NA_real_, vcount(g))
)

power_result <- tryCatch(
  power_centrality(
    g,
    exponent = 0.01,
    weights = NA
  ),
  error = function(e) rep(NA_real_, vcount(g))
)

subgraph_result <- tryCatch(
  subgraph_centrality(g),
  error = function(e) rep(NA_real_, vcount(g))
)

centrality_table$alpha <- alpha_result
centrality_table$power <- power_result
centrality_table$subgraph <- subgraph_result

# HITS: hub + authority.
hits_result <- hits_scores(
  g,
  scale = TRUE
)

centrality_table$hub <- hits_result$hub
centrality_table$authority <- hits_result$authority


unit_weights <- rep(1, ecount(g))

centrality_table$strength <- strength(
  g,
  weights = unit_weights
)

cat("\nЦентральности вершин:\n")
print(centrality_table)

cat("\nЦентральность рёбер: edge betweenness\n")

edge_betweenness_result <- edge_betweenness(
  g,
  directed = FALSE,
  weights = NA
)

edge_centrality_table <- data.frame(
  edge_id = seq_len(ecount(g)),
  from = ends(g, E(g), names = FALSE)[, 1],
  to = ends(g, E(g), names = FALSE)[, 2],
  edge_betweenness = edge_betweenness_result
)

print(edge_centrality_table)

cat("\nДругие центральности:\n")
cat("- diversity рассчитывается отдельно ниже, поскольку для неё нужны веса рёбер.\n")


# Diversity требует весов рёбер.
set.seed(42)
diversity_weights <- runif(
  ecount(g),
  min = 0.5,
  max = 1.5
)

diversity_result <- diversity(
  g,
  weights = diversity_weights
)

cat("\nDiversity центральность вершин:\n")
print(diversity_result)


# ============================================================
# 10. Кластеры
# ============================================================

cat("\n========== 8. СООБЩЕСТВА / КЛАСТЕРЫ ==========\n")

# Для каждого алгоритма выводится:
# - количество кластеров;
# - принадлежность вершин;
# - модульность, если она рассчитана.

print_community <- function(name, community_object) {

  cat("\n---", name, "---\n")

  print(community_object)

  cat("Количество кластеров:", length(sizes(community_object)), "\n")

  cat("Размеры кластеров:\n")
  print(sizes(community_object))

  cat("Membership:\n")
  print(membership(community_object))

  cat("Modularity:\n")

  mod_value <- tryCatch(
    modularity(community_object),
    error = function(e) NA_real_
  )

  print(mod_value)
}


# 1. Edge betweenness
community_edge_betweenness <- cluster_edge_betweenness(g)
print_community(
  "Edge betweenness",
  community_edge_betweenness
)


# 2. Fast greedy
community_fast_greedy <- cluster_fast_greedy(g)
print_community(
  "Fast greedy",
  community_fast_greedy
)


# 3. Fluid communities
community_fluid <- cluster_fluid_communities(
  g,
  no = 3
)
print_community(
  "Fluid communities",
  community_fluid
)


# 4. Infomap
community_infomap <- cluster_infomap(g)
print_community(
  "Infomap",
  community_infomap
)


# 5. Label propagation
community_label_prop <- cluster_label_prop(g)
print_community(
  "Label propagation",
  community_label_prop
)


# 6. Leading eigenvector
community_leading_eigen <- cluster_leading_eigen(g)
print_community(
  "Leading eigenvector",
  community_leading_eigen
)


# 7. Leiden
community_leiden <- cluster_leiden(
  g,
  objective_function = "modularity"
)
print_community(
  "Leiden",
  community_leiden
)


# 8. Louvain
community_louvain <- cluster_louvain(g)
print_community(
  "Louvain",
  community_louvain
)


# 9. Spinglass
community_spinglass <- cluster_spinglass(g)
print_community(
  "Spinglass",
  community_spinglass
)


# 10. Walktrap
community_walktrap <- cluster_walktrap(g)
print_community(
  "Walktrap",
  community_walktrap
)


# 11. Optimal community structure.
# Точный поиск может работать существенно дольше остальных
# методов, поэтому его можно отключить переключателем.
RUN_CLUSTER_OPTIMAL <- TRUE

if (RUN_CLUSTER_OPTIMAL) {

  community_optimal <- tryCatch(
    cluster_optimal(g),
    error = function(e) {
      message(
        "cluster_optimal() не выполнился: ",
        conditionMessage(e)
      )
      NULL
    }
  )

  if (!is.null(community_optimal)) {

    print_community(
      "Optimal",
      community_optimal
    )
  }
}


# ============================================================
# 11. Ядра графа
# ============================================================

cat("\n========== 9. ЯДРА ГРАФА ==========\n")

# Coreness каждой вершины.
coreness_result <- coreness(
  g,
  mode = "all"
)

cat("Coreness вершин:\n")
print(coreness_result)

# Максимальное k.
max_k <- max(coreness_result)

cat("\nМаксимальное значение k-core:", max_k, "\n")

# Сами подграфы k-core.
k_core_graphs <- list()

for (k in 1:max_k) {

  core_vertices <- V(g)[coreness_result >= k]

  k_core_graphs[[paste0("k_", k)]] <- induced_subgraph(
    g,
    core_vertices
  )

  cat(
    "\nk =", k,
    ": вершин =", vcount(k_core_graphs[[paste0("k_", k)]]),
    ", рёбер =", ecount(k_core_graphs[[paste0("k_", k)]])
  )
  cat("\n")
}


# ============================================================
# 12. Все простые циклы
# ============================================================

cat("\n========== 10. ПРОСТЫЕ ЦИКЛЫ ==========\n")

cycles_3_6 <- simple_cycles(
  g,
  min = 3,
  max = 6
)

cat(
  "Количество найденных простых циклов длиной 3-6:",
  length(cycles_3_6),
  "\n"
)

# Вывести циклы в удобном виде.
if (length(cycles_3_6) > 0) {

  for (i in seq_along(cycles_3_6)) {

    cat(
      i, ": ",
      paste(
        as_ids(cycles_3_6[[i]]),
        collapse = " -> "
      ),
      "\n"
    )
  }
}


# ============================================================
# 13. Матрица смежности
# ============================================================

cat("\n========== 11. МАТРИЦЫ ==========\n")

adjacency_matrix <- as_adjacency_matrix(
  g,
  sparse = FALSE
)

cat("\nМатрица смежности:\n")
print(adjacency_matrix)


# ============================================================
# 14. Матрица инцидентности
# ============================================================

incidence_matrix <- matrix(
  0,
  nrow = vcount(g),
  ncol = ecount(g)
)

edge_ends <- ends(
  g,
  E(g),
  names = FALSE
)

for (edge_id in seq_len(ecount(g))) {

  v1 <- edge_ends[edge_id, 1]
  v2 <- edge_ends[edge_id, 2]

  incidence_matrix[v1, edge_id] <- 1
  incidence_matrix[v2, edge_id] <- 1
}

rownames(incidence_matrix) <- V(g)$name

colnames(incidence_matrix) <- paste0(
  "e",
  seq_len(ecount(g))
)

cat("\nМатрица инцидентности:\n")
print(incidence_matrix)


# ============================================================
# 15. Матрица модулярности
# ============================================================

# В качестве разбиения используем Leiden.
membership_leiden <- membership(
  community_leiden
)

modularity_matrix_result <- modularity_matrix(
  g,
  membership = membership_leiden,
  weights = NULL,
  resolution = 1,
  directed = FALSE
)

cat("\nМатрица модулярности:\n")
print(modularity_matrix_result)

cat(
  "\nМодульность разбиения Leiden:",
  modularity(community_leiden),
  "\n"
)


# ============================================================
# 16. Эксцентриситеты вершин
# ============================================================

cat("\n========== 12. ЭКСЦЕНТРИСИТЕТЫ ==========\n")

eccentricity_result <- eccentricity(
  g,
  mode = "all",
  weights = NA
)

eccentricity_table <- data.frame(
  vertex = seq_len(vcount(g)),
  name = V(g)$name,
  eccentricity = eccentricity_result
)

print(eccentricity_table)


# ============================================================
# 17. Радиус, диаметр и центр
# ============================================================

cat("\n========== 13. РАДИУС, ДИАМЕТР И ЦЕНТР ==========\n")

graph_radius <- radius(
  g,
  mode = "all",
  weights = NA
)

graph_diameter <- diameter(
  g,
  directed = FALSE,
  weights = NA
)

graph_center <- which(
  eccentricity_result == graph_radius
)

cat("Радиус графа:", graph_radius, "\n")
cat("Диаметр графа:", graph_diameter, "\n")

cat("Центр графа (ID вершин):\n")
print(graph_center)

cat("Центр графа (имена вершин):\n")
print(V(g)$name[graph_center])