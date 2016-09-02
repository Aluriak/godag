OUTPUT=output.json
IDNAME=mapping.json
USE_ID=--use_id

all:
	python godag.py $(OUTPUT) $(IDNAME) $(USE_ID)
