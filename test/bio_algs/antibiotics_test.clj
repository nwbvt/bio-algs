(ns bio-algs.antibiotics-test
  (:require [clojure.test :refer :all]
            [bio-algs.antibiotics :refer :all]))

(deftest protein-translation
  (testing "That dna translates into proteins correctly"
    (is (= "MAMAPRTEINSTRING"
           (translate "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA")))))

(deftest peptide-encoding
  (testing "Finding encodings for a peptide"
    (is (= ["ATGGCC" "GGCCAT" "ATGGCC"]
           (find-encoding "ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA" "MA")))))
