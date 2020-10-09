/****************************************************************************************** 
 *	Chili DirectX Framework Version 16.07.20											  *	
 *	Game.h																				  *
 *	Copyright 2016 PlanetChili.net <http://www.planetchili.net>							  *
 *																						  *
 *	This file is part of The Chili DirectX Framework.									  *
 *																						  *
 *	The Chili DirectX Framework is free software: you can redistribute it and/or modify	  *
 *	it under the terms of the GNU General Public License as published by				  *
 *	the Free Software Foundation, either version 3 of the License, or					  *
 *	(at your option) any later version.													  *
 *																						  *
 *	The Chili DirectX Framework is distributed in the hope that it will be useful,		  *
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of						  *
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the						  *
 *	GNU General Public License for more details.										  *
 *																						  *
 *	You should have received a copy of the GNU General Public License					  *
 *	along with The Chili DirectX Framework.  If not, see <http://www.gnu.org/licenses/>.  *
 ******************************************************************************************/
#pragma once

#include "Keyboard.h"
#include "Mouse.h"
#include "Graphics.h"
#include "FrameTimer.h"

class Game
{
public:
	Game( class MainWindow& wnd );
	Game( const Game& ) = delete;
	Game& operator=( const Game& ) = delete;
	~Game();
	void Go();
private:
	void ComposeFrame();
	void UpdateModel();
	/********************************/
	/*  User Functions              */
	int GetId(int i, int j);
	void DrawDensity();
	void DrawVelocities(bool separated);
	void DensitySolver(float brushAmountPerSec, float brushRadius, float diffusionRate, float dt);
	void AddDensity(float AmountPerSec, float radius, float dt);
	void Diffusion(float diffusionRate, float dt);
	/********************************/
private:
	MainWindow& wnd;
	Graphics gfx;
	/********************************/
	/*  User Variables              */
	FrameTimer ft;
	static constexpr int n = 41;
	static constexpr int N = n + 2;
	static constexpr float cellDimension = float(Graphics::ScreenWidth) / float(N);
	float* density = nullptr;
	float* prev_density = nullptr;
	Vec2* velocity = nullptr;
	Vec2* prev_velocity = nullptr;
	static constexpr float diffusionRate = 1.0f;
	static constexpr float DPS = 1.0f;
	static constexpr float brushRadius = 1.0f;
	//test
	/********************************/
};