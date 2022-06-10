#For source code instalations do
KAKA = /home/live/repos/MFEM/mfem
CONFIG_MK = $(KAKA)/build/config/config.mk
GENERAL = $(KAKA)/general

#For Spack instalations do
#MFEM_INSTALL_DIR = /opt/spack/opt/spack/linux-pop21-icelake/gcc-11.2.0/mfem-4.3.0-wor2wndreilv3jrwzqxm5a54cvqgrvaq
#CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk
#GENERAL = $(MFEM_INSTALL_DIR)/include/mfem/general

SHARE_DIR = NULL
PROCCESORS = 1

include $(CONFIG_MK)
