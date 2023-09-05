module resource_module
  !Module containing the data type for resouce management.

#ifdef _OPENACC
  use openacc
#endif
  
  implicit none
  
  public resource
  
  type resource
     integer :: num_cpus
     !! Number of bare cpus (not a gpu manager)
     integer :: num_gpus
     !! Number of gpus
     integer :: this_node
     !! ID of node where this image resides
     logical :: gpu_manager
     !! Is this image a gpu manager?
     character*(30) :: gpu_name = trim("None"), &
          gpu_vendor = trim("None"), gpu_driver = trim("None")
     !! Information about gpu device
     
   contains
     procedure :: initialize, report
  end type resource

contains

  subroutine initialize(self)
    class(resource), intent(out) :: self

    integer :: igpus, im, iset
    integer :: devicenum
    integer(acc_device_property):: property

    character*(100) :: string
    character(len=10), allocatable :: hostname(:)[:]
    character(len=10), allocatable :: hostname_set(:)

    self%num_gpus = 0
#ifdef _OPENACC
    self%num_gpus = acc_get_num_devices(acc_device_default)
#endif
    self%num_cpus = num_images() - self%num_gpus
    
    allocate(hostname(num_images())[*])
    call hostnm(hostname(this_image()))
    if(this_image() == 1) then
       do im = 1, num_images()
          hostname(im) = hostname(im)[im]
       end do
       print*, 'hostname: ', hostname(:)
    end if
    sync all
    call co_broadcast(hostname, 1)
    sync all
    
    call create_set(hostname, hostname_set)

    self%this_node = findloc(hostname_set, hostname(this_image()), 1)
    self%gpu_manager = &
         (this_image() == findloc(hostname, hostname_set(self%this_node), 1)) &
         .and. self%num_gpus > 0

#ifdef _OPENACC
    if(self%gpu_manager) then
       do igpus = 0, self%num_gpus - 1 !Mind the 0 based indexing of openacc
          property = acc_property_name
          call acc_get_property_string(igpus, acc_get_device_type(), &
               property, string)
          self%gpu_name = trim(string)

          property = acc_property_vendor
          call acc_get_property_string(igpus, acc_get_device_type(), &
               property, string)
          self%gpu_vendor = trim(string)

          property = acc_property_driver
          call acc_get_property_string(igpus, acc_get_device_type(), &
               property, string)
          self%gpu_driver = trim(string)
       end do
    end if
#endif
  end subroutine initialize

  subroutine report(self)
    class(resource), intent(in) :: self

    if(self%gpu_manager) then
       print*, self%this_node, this_image(), self%gpu_name, &
            self%gpu_vendor, self%gpu_driver
    end if
  end subroutine report
end module resource_module
