(ns bio-algs.proteomics-test  
  (:require [clojure.test :refer :all]
            [bio-algs.proteomics :refer :all]
            [bio-algs.core :refer [parse-graph]]))

(deftest ideal-spectrums
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
                          "429->486:G"] :edge-parser identity))))
  (testing "Generating a spectrum from a peptide"
    (is (= (sort (gen-spectrum "ANFPG")) [0 57 71 154 185 301 332 415 429 486])))
  (testing "the explains function works"
    (is (true? (explains-spectrum? [57 71 154 185 301 332 415 429 486] "GPFNA")))
    (is (false? (explains-spectrum? [57 71 154 185 301 332 415 429] "GPFNA")))
    (is (false? (explains-spectrum? [57 71 154 185 301 332 415 429 486 512] "GPFNA"))))
  (testing "Using the graph to decode the spectrum"
    (is (= (decode-spec-graph [57 71 154 185 301 332 415 429 486])
           "ANFPG"))))

(deftest noisy-spectrums
  (testing "Converting a peptide to a peptide vector and back again"
    (let [peptide "TFPRGPHSPRVVDIRCCQQMNDHQSIDWQYSIYFM"   
          v (peptide-vector peptide)]
      (is (= (count v) 4234))
      (is (= (reduce + v) 35))
      (is (= (subvec (vec v) 80 120)
             [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]))
      (is (= (from-pv v) peptide)))))
