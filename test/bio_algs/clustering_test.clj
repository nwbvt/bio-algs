(ns bio-algs.clustering-test
  (:require [clojure.test :refer :all]
            [bio-algs.clustering :refer :all]
            [clojure.math.numeric-tower :refer [abs round expt]]))

(defn round-all
  "round the numbers in the given data structure m to n decimal points"
  [n m]
  (if (seq? m) (map (partial round-all n) m)
    (if (number? m)
      (let [d (expt 10 n)]
        (double (/ (round (* m d)) d)))
      m)))

(deftest farthest-first-clustering
  (testing "the farthest first clustering algorithm"
    (is (= (ff-clusters 3 2 [[0.0 0.0]
                             [5.0 5.0]
                             [0.0 5.0]
                             [1.0 1.0]
                             [2.0 2.0]
                             [3.0 3.0]
                             [1.0 2.0]])
           [[0.0 0.0]
            [5.0 5.0]
            [0.0 5.0]]))))

(deftest k-means-clustering
  (testing "computing the squared error distortion"
    (is (< (abs (- (squared-error 2 [[2.31 4.55] [5.96 9.08]]
                                  [[3.42 6.03]
                                   [6.23 8.25]
                                   [4.76 1.64]
                                   [4.47 4.33]
                                   [3.95 7.61]
                                   [8.93 2.97]
                                   [9.74 4.03]
                                   [1.73 1.28]
                                   [9.72 5.01]
                                   [7.27 3.77]])
                   18.246))
           0.001)))
  (testing "The Lloyd k-means algorithm"
    (is (= (round-all 3 (k-means 2 2 [[1.3 1.1]
                                      [1.3 0.2]
                                      [0.6 2.8]
                                      [3.0 3.2]
                                      [1.2 0.7]
                                      [1.4 1.6]
                                      [1.2 1.0]
                                      [1.2 1.1]
                                      [0.6 1.5]
                                      [1.8 2.6]
                                      [1.2 1.3]
                                      [1.2 1.0]
                                      [0.0 1.9]]))
           [[1.06 1.14]
            [1.8 2.867]]))))

(deftest soft-clustering
  (testing "clustering using expectation maximization"
    (is (= (round-all 3 (e-m-cluster 2 2 2.7
                                     [[1.3 1.1]
                                      [1.3 0.2]
                                      [0.6 2.8]
                                      [3.0 3.2]
                                      [1.2 0.7]
                                      [1.4 1.6]
                                      [1.2 1.0]
                                      [1.2 1.1]
                                      [0.6 1.5]
                                      [1.8 2.6]
                                      [1.2 1.3]
                                      [1.2 1.0]
                                      [0.0 1.9]] 100))
           [[1.662 2.623]
            [1.075 1.148]]))))

(deftest hierarchical-clustering
  (testing "clustering using hierarchical clustering"
    (is (= (h-clusters 7
             [[0.00 0.74 0.85 0.54 0.83 0.92 0.89]
              [0.74 0.00 1.59 1.35 1.20 1.48 1.55]
              [0.85 1.59 0.00 0.63 1.13 0.69 0.73]
              [0.54 1.35 0.63 0.00 0.66 0.43 0.88]
              [0.83 1.20 1.13 0.66 0.00 0.72 0.55]
              [0.92 1.48 0.69 0.43 0.72 0.00 0.80]
              [0.89 1.55 0.73 0.88 0.55 0.80 0.00]])
           [[4 6]
            [5 7]
            [3 4 6]
            [1 2]
            [3 4 5 6 7]
            [1 2 3 4 5 6 7]]))))
