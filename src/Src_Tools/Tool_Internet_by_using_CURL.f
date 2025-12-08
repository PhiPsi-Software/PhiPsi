 
      subroutine Tool_Internet_by_using_CURL
      
      use Global_Float_Type
      use Global_Common
      
      implicit none
      character internetget*200,url*80,outputfile*80
      
      url = 'phipsi.top/run_phipsi_visit.html'
      
      outputfile  = trim(PhiPsi_Current_Directory) 
     &       //'\PhiPsi_visit_phipsi_website.log'
      internetget = trim('curl '//url// '-o '//outputfile)
      call system(internetget)
      
      
      return 
      end SUBROUTINE Tool_Internet_by_using_CURL                        
