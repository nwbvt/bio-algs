(ns bio-algs.clustering-test
  (:require [clojure.test :refer :all]
            [bio-algs.clustering :refer :all]
            [clojure.math.numeric-tower :refer [abs round]]))

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
    (let [results (k-means 2 2 [[1.3 1.1]
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
                                [0.0 1.9]])]
      (is (= (for [center results] (for [i center] (/ (round (* 1000 i)) 1000.0)))
             [[1.060 1.140]
              [1.800 2.867]])))))
