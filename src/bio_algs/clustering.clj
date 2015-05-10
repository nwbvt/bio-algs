(ns bio-algs.clustering
  (:require [bio-algs.core :refer [read-file to-floats to-ints]]
            [clojure.string :refer [join]]
            [clojure.math.numeric-tower :refer [expt sqrt]]))

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
