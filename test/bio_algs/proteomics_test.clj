(ns bio-algs.proteomics-test  
  (:require [clojure.test :refer :all]
            [bio-algs.proteomics :refer :all]
            [bio-algs.core :refer [parse-graph]]))

(deftest spectrums
  (testing "Find the graph for a spectrum"
    (is (= (spec-graph [57 71 154 185 301 332 415 429 486])
           (parse-graph [ "0->57:G"
                          "0->71:A"
                          "57->154:P"
                          "57->185:Q"
                          "71->185:N"
                          "154->301:F"
                          "185->332:F"
                          "301->415:N"
                          "301->429:Q"
                          "332->429:P"
                          "415->486:A"
                          "429->486:G"] :edge-parser identity)))))
