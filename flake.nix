{ 
  description = "Julia env flake";

  inputs = {
    nixgl.url = "github:nix-community/nixGL";
    nixpkgs.url = "github:nixos/nixpkgs?ref=nixos-unstable";
    scientific-fhs.url =  "github:olynch/scientific-fhs";
  };
  outputs = { self, nixgl, nixpkgs, scientific-fhs }:
    let 
      system = "x86_64-linux";
      pkgs = import nixpkgs {
        inherit system; 
        overlays = [nixgl.overlay];
      };
    in {
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
