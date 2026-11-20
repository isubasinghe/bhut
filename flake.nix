{
  description = "Barnes-Hut N-body simulation with MPI, and a video of the result";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

  outputs = { self, nixpkgs }:
    let
      systems = [ "x86_64-linux" "aarch64-linux" "x86_64-darwin" "aarch64-darwin" ];
      forAllSystems = f: nixpkgs.lib.genAttrs systems (system: f nixpkgs.legacyPackages.${system});
    in
    {
      # Enter with `nix develop`, then run `just`.
      devShells = forAllSystems (pkgs: {
        default = pkgs.mkShell {
          packages = [
            pkgs.openmpi # mpicxx, mpirun
            pkgs.just
            pkgs.ffmpeg
            pkgs.rustc
            pkgs.cargo
            pkgs.vulkan-loader
            (pkgs.python3.withPackages (ps: [ ps.numpy ps.matplotlib ]))
          ];
          shellHook = ''
            export LD_LIBRARY_PATH="${pkgs.lib.makeLibraryPath [ pkgs.vulkan-loader ]}''${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
          '';
        };
      } // pkgs.lib.optionalAttrs pkgs.stdenv.hostPlatform.isLinux (
        let
          cudaPkgs = import nixpkgs {
            system = pkgs.stdenv.hostPlatform.system;
            config = { allowUnfree = true; cudaCapabilities = [ "7.5" ]; };
          };
          cuda = cudaPkgs.cudaPackages_12_9;
        in {
          cuda = (cudaPkgs.mkShell.override { stdenv = cuda.backendStdenv; }) {
            packages = [ cuda.cuda_nvcc cuda.cuda_cudart cuda.cccl cudaPkgs.just cudaPkgs.python3 ];
            # CUDA components have separate Nix outputs; give nvcc explicit paths.
            NVCC_PREPEND_FLAGS = builtins.concatStringsSep " " [
              "-ccbin ${cuda.backendStdenv.cc}/bin/c++"
              "-I${cudaPkgs.lib.getDev cuda.cuda_cudart}/include"
              "-I${cudaPkgs.lib.getDev cuda.cccl}/include"
              "-L${cudaPkgs.lib.getLib cuda.cuda_cudart}/lib"
              "--cudart=shared"
              "--cudadevrt=none"
            ];
            shellHook = ''
              export LD_LIBRARY_PATH="/run/opengl-driver/lib:${cudaPkgs.lib.makeLibraryPath [ cuda.cuda_cudart ]}''${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
            '';
          };
        }
      ));
    };
}
