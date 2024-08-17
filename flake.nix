{ 
  description = "Julia ScientificFHS flake";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs?ref=nixos-unstable";
    scientific-fhs.url =  "github:olynch/scientific-fhs";
  };
  outputs = { self, nixpkgs, scientific-fhs }:
    let 
      system = "x86_64-linux";
      pkgs = import nixpkgs {
        inherit system; 
        overlays = [];
      };
    in {
      imports = [scientific-fhs.nixosModules.default];
      programs.scientific-fhs = {
        enable = true;
        juliaVersions = [
          {
            version = "julia_18";
            default = true;
          }
        ];
        enableNVIDIA = false;
      };
      devShells.${system}.default = pkgs.mkShell {
        packages = with pkgs; [
          julia-bin
          gr-framework
          stdenv.cc.cc.lib qt5.qtbase qt5Full libGL
          glxinfo
          glfw
        ];

      };
      
    };
}
