OUTPUT=output/godag
IDNAME=output/mapping
USE_ID=--use_id

json:
	python godag.py $(OUTPUT).json $(IDNAME).json $(USE_ID) --json

dsv:
	python godag.py $(OUTPUT).dsv $(IDNAME).dsv $(USE_ID)
