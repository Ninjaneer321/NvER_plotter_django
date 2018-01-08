from django import forms

class Gene_searchForm(forms.Form):
    gene_name = forms.CharField(max_length=20)

class NvERTxForm(forms.Form):
    nvertx_1 = forms.CharField(max_length=25)
    nvertx_2 = forms.CharField(max_length=25, required=False)
    nvertx_3 = forms.CharField(max_length=25, required=False)
    nvertx_4 = forms.CharField(max_length=25, required=False)
    nvertx_5 = forms.CharField(max_length=25, required=False)
    log2 = forms.BooleanField(initial=False, required=False)

class ConvertForm(forms.Form):
	nvertx = forms.CharField(max_length=25)

class comparisonForm(forms.Form):
	CHOICES= (
	('regen_0_2','Regeneration: 0 to 2 Hpa'),
	('regen_0_4','Regeneration: 0 to 4 Hpa'),
	('regen_0_8','Regeneration: 0 to 8 Hpa'),
	('regen_0_12','Regeneration: 0 to 12 Hpa'),
	('regen_0_16','Regeneration: 0 to 16 Hpa'),
	('regen_0_20','Regeneration: 0 to 20 Hpa'),
	('regen_0_24','Regeneration: 0 to 24 Hpa'),
	('regen_0_36','Regeneration: 0 to 36 Hpa'),
	('regen_0_48','Regeneration: 0 to 48 Hpa'),
	('regen_0_60','Regeneration: 0 to 60 Hpa'),
	('regen_0_72','Regeneration: 0 to 72 Hpa'),
	('regen_0_96','Regeneration: 0 to 96 Hpa'),
	('regen_0_120','Regeneration: 0 to 120 Hpa'),
	('regen_0_144','Regeneration: 0 to 144 Hpa'),
	('fischer_7_8','Embryogenesis (Fischer): 7 to 8 Hpf'),
	('fischer_7_9','Embryogenesis (Fischer): 7 to 9 Hpf'),
	('fischer_7_10','Embryogenesis (Fischer): 7 to 10 Hpf'),
	('fischer_7_11','Embryogenesis (Fischer): 7 to 11 Hpf'),
	('fischer_7_12','Embryogenesis (Fischer): 7 to 12 Hpf'),
	('fischer_7_13','Embryogenesis (Fischer): 7 to 13 Hpf'),
	('fischer_7_14','Embryogenesis (Fischer): 7 to 14 Hpf'),
	('fischer_7_15','Embryogenesis (Fischer): 7 to 15 Hpf'),
	('fischer_7_16','Embryogenesis (Fischer): 7 to 16 Hpf'),
	('fischer_7_17','Embryogenesis (Fischer): 7 to 17 Hpf'),
	('fischer_7_18','Embryogenesis (Fischer): 7 to 18 Hpf'),
	('fischer_7_19','Embryogenesis (Fischer): 7 to 19 Hpf'),
	('helm_7_12','Embryogenesis (Helm): 7 to 12 Hpf'),
	('helm_7_24','Embryogenesis (Helm): 7 to 24 Hpf'),
	('helm_7_120','Embryogenesis (Helm): 7 to 120 Hpf'),
	('helm_7_240','Embryogenesis (Helm): 7 to 240 Hpf'),
	('warner_24_48','Embryogenesis (Warner): 24 to 48 Hpf'),
	('warner_24_72','Embryogenesis (Warner): 24 to 72 Hpf'),
	('warner_24_96','Embryogenesis (Warner): 24 to 96 Hpf'),
	('warner_24_120','Embryogenesis (Warner): 24 to 120 Hpf'),
	('warner_24_144','Embryogenesis (Warner): 24 to 144 Hpf'),
	('warner_24_168','Embryogenesis (Warner): 24 to 168 Hpf'),
	('warner_24_192','Embryogenesis (Warner): 24 to 192 Hpf'),
	)
	comparison = forms.ChoiceField(widget=forms.Select, choices=CHOICES)