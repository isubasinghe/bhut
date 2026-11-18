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
            (pkgs.python3.withPackages (ps: [ ps.numpy ps.matplotlib ]))
          ];
        };
      });
    };
}
