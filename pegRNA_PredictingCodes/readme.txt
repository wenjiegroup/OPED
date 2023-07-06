this is for prediction only, 'predict_efficiency_of_ClinVar.py' is enough.

run in shell:
	python predict_efficiency_of_ClinVar.py

default input file: Data/pegRNA.GRCh38.Insertion.txt


actually, "Target(47bp)	PBS	RT" are used, users should provide Target \ PBS \ RT


Trained model "Model_Trained/pegRNA_Model_Merged_saved.order3_decoder.pt" are used.


Output stored in "Output/Efficiency.pegRNA.GRCh38.Insertion.txt", that is to be used for users.
