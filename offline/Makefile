ifndef GO4SYS
GO4SYS = $(shell go4-config --go4sys)
endif

include $(GO4SYS)/Makefile.config

## normally should be like this for every module, but can be specific

ifdef GO4PACKAGE
NFISPTB_DIR         = Go4nfisPTB
else
NFISPTB_DIR         = .
endif

NFISPTB_LINKDEF     = $(NFISPTB_DIR)/nfisLinkDef.$(HedSuf)
NFISPTB_LIBNAME     = $(GO4_USERLIBNAME)

## must be similar for every module

NFISPTB_DICT        = $(NFISPTB_DIR)/$(DICT_PREFIX)nfis
NFISPTB_DH          = $(NFISPTB_DICT).$(HedSuf)
NFISPTB_DS          = $(NFISPTB_DICT).$(SrcSuf)
NFISPTB_DO          = $(NFISPTB_DICT).$(ObjSuf)

NFISPTB_H           = $(filter-out $(NFISPTB_DH) $(NFISPTB_LINKDEF), $(wildcard $(NFISPTB_DIR)/*.$(HedSuf)))
NFISPTB_S           = $(filter-out $(NFISPTB_DS), $(wildcard $(NFISPTB_DIR)/*.$(SrcSuf)))
NFISPTB_O           = $(NFISPTB_S:.$(SrcSuf)=.$(ObjSuf))

NFISPTB_DEP         =  $(NFISPTB_O:.$(ObjSuf)=.$(DepSuf))
NFISPTB_DDEP        =  $(NFISPTB_DO:.$(ObjSuf)=.$(DepSuf))

NFISPTB_LIB         =  $(NFISPTB_DIR)/$(NFISPTB_LIBNAME).$(DllSuf)

# used in the main Makefile

EXAMPDEPENDENCS    += $(NFISPTB_DEP) $(NFISPTB_DDEP)

ifdef DOPACKAGE
DISTRFILES         += $(NFISPTB_S) $(NFISPTB_H) $(NFISPTB_LINKDEF)
DISTRFILES         += $(NFISPTB_DIR)/Readme.txt
DISTRFILES         += $(wildcard $(NFISPTB_DIR)/*.C)
endif

##### local rules #####

all::     $(NFISPTB_LIB)

$(NFISPTB_LIB):  $(NFISPTB_O) $(NFISPTB_DO) $(ANAL_LIB_DEP)
	@$(MakeLibrary) $(NFISPTB_LIBNAME) "$(NFISPTB_O) $(NFISPTB_DO)" $(NFISPTB_DIR) $(NFISPTB_LINKDEF) "$(ANAL_LIB_DEP)" $(NFISPTB_DS) "$(NFISPTB_H)"

$(NFISPTB_DS): $(NFISPTB_H)  $(NFISPTB_LINKDEF)
	@$(ROOTCINTGO4) $(NFISPTB_LIB) $(NFISPTB_H) $(NFISPTB_LINKDEF)

clean-bin::
	@rm -f $(NFISPTB_O) $(NFISPTB_DO)
	@rm -f $(NFISPTB_DEP) $(NFISPTB_DDEP) $(NFISPTB_DS) $(NFISPTB_DH)

clean:: clean-bin
	@$(CleanLib) $(NFISPTB_LIBNAME) $(NFISPTB_DIR)

include $(GO4SYS)/Makefile.rules
