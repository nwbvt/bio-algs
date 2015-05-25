(ns bio-algs.clustering
  (:require [bio-algs.core :refer [read-file to-floats to-ints]]
            [clojure.string :refer [join]]
            [clojure.math.numeric-tower :refer [expt sqrt]]
            [clojure.set :refer [union]]))

(defn d
  "the minimum distance from m dimensional point x to the list of points"
  [m points x]
  (apply min (for [y points]
               (sqrt (apply + (for [i (range m) :let [xi (nth x i 0) yi (nth y i 0)]]
                                (expt (- xi yi) 2)))))))

(defn ff-clusters
  "use furthest first clustering to find k clusters for the given points in m dimensional space"
  [k m points]
  (loop [clusters [(first points)]]
    (if (= k (count clusters))
      clusters
      (recur (conj clusters (apply max-key (partial d m clusters) points))))))

(defn run-ff-clustering
  [input-file]
  (let [[[k m] & data] (to-floats (read-file input-file))
        clusters (ff-clusters (int k) (int m) data)]
    (map (partial join " ") clusters)))

(defn squared-error
  "Find the squared error of the clusters for the data in m dimensional space"
  [m clusters data]
  (/ (apply + (for [x data] (expt (d m clusters x) 2))) (count data)))

(defn run-squared-error
  [input-file]
  (let [input (read-file input-file)
        [[k m]] (to-ints [(first input)])
        clusters (to-floats (take k (rest input)))
        data (to-floats (drop (+ 2 k) input))]
    (squared-error m clusters data)))
(ns bio-algs.clustering
  (:require [bio-algs.core :refer [read-file to-floats]]
            [clojure.string :refer [join]]
            [clojure.math.numeric-tower :refer [expt sqrt]]
            ))

(defn d
  "the minimum distance from m dimensional point x to the list of points"
  [m points x]
  (apply min (for [y points]
               (sqrt (apply + (for [i (range m) :let [xi (nth x i 0) yi (nth y i 0)]]
                                (expt (- xi yi) 2)))))))

(defn ff-clusters
  "use furthest first clustering to find k clusters for the given points in m dimensional space"
  [k m points]
  (loop [clusters [(first points)]]
    (if (= k (count clusters))
      clusters
      (recur (conj clusters (apply max-key (partial d m clusters) points))))))

(defn run-ff-clustering
  [input-file]
  (let [[[k m] & data] (to-floats (read-file input-file))
        clusters (ff-clusters (int k) (int m) data)]
    (map (partial join " ") clusters)))

(defn squared-error
  "Find the squared error of the clusters for the data in m dimensional space"
  [m clusters data]
  (/ (apply + (for [x data] (expt (d m clusters x) 2))) (count data)))

(defn run-squared-error
  [input-file]
  (let [input (read-file input-file)
        [[k m]] (to-ints [(first input)])
        clusters (to-floats (take k (rest input)))
        data (to-floats (drop (+ 2 k) input))]
    (squared-error m clusters data)))

(defn- find-best-cluster
  [m centers x]
  (apply min-key #(d m [%] x) centers))

(defn center-of-gravity
  [m points]
  (for [i (range m)]
    (/ (apply + (map #(nth % i) points)) (count points))))

(defn k-means
  "Run lloyd's k-means algorithm to find the centers of the clusters"
  [k m points]
  (loop [centers (take k points)]
    (let [clusters (group-by (partial find-best-cluster m centers) points)
          new-centers (map (partial center-of-gravity m) (vals clusters))]
      (if (= (set new-centers) (set centers)) centers
        (recur new-centers)))))

(defn run-k-means
  [input-file]
  (let [[[k m] & data] (to-floats (read-file input-file))
        centers (k-means (int k) (int m) data)]
    (map (partial join " ") centers)))

(defn calc-responsibilities
 [m beta centers]
 (fn [point]
   (let [ds (map #(d m [%] point) centers)
         numerators (map #(expt Math/E (* (- beta) %)) ds)
         total (apply + numerators)]
     (map #(/ % total) numerators))))

(defn find-center
  [data m responsibilities]
  (let [total (apply + responsibilities)]
    (for [i (range m)] (/ (apply + (pmap #(* (nth responsibilities %) (nth (nth data %) i)) (range (count data)))) total))))

(defn e-m-cluster
  [k m beta data n]
  (loop [centers (take k data) remaining n]
    (if (zero? remaining) centers
      (let [soft-clusters (pmap (calc-responsibilities m beta centers) data)
            new-centers (for [i (range k)] (find-center data m (map #(nth % i) soft-clusters)))]
        (recur new-centers (dec remaining))))))

(defn run-e-m-cluster
  [input-file]
  (let [[[k m] [beta] & data] (to-floats (read-file input-file))
        centers (e-m-cluster (int k) (int m) beta data 100)]
    (map (partial join " ") centers)))

(defn cluster-dist
  "Return the average distance between the cluster members"
  [dist-matrix [cluster1 cluster2]]
  (/ (apply + (for [p1 cluster1 p2 cluster2] (-> dist-matrix (nth (dec p1)) (nth (dec p2))))) (* (count cluster1) (count cluster2))))

(defn closest-clusters
  "Finds the pair of clusters that are the closest"
  [dist-matrix clusters]
  (apply min-key (partial cluster-dist dist-matrix) (for [cluster1 clusters cluster2 clusters :when (not= cluster1 cluster2)] [cluster1 cluster2])))

(defn next-clusters
  [dist-matrix clusters]
  (let [[cluster1 cluster2] (closest-clusters dist-matrix clusters)]
    (conj (remove #(or (= cluster1 %) (= cluster2 %)) clusters) (union cluster1 cluster2))))

(defn gen-h-clusters
  "Use the hierarchical clustering algorithm to generate a sequence of clusters
   Where the newest cluster is always the first member"
  [n dist-matrix]
  (let [initial-clusters (map #(set [%]) (range 1 (inc n)))]
    (take n (iterate (partial next-clusters dist-matrix) initial-clusters))))

(defn h-clusters
  "List the clusters generated by the hierarchical clustering algorithm"
  [n dist-matrix]
  (map #(sort (seq %)) (map first (rest (gen-h-clusters n dist-matrix)))))

(defn run-h-clusters
  [input-file]
  (let [[[n] & dist-matrix] (to-floats (read-file input-file))
        clusters (h-clusters (int n) dist-matrix)]
    (map (partial join " ") clusters)))
