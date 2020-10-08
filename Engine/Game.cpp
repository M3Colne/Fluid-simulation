/****************************************************************************************** 
 *	Chili DirectX Framework Version 16.07.20											  *	
 *	Game.cpp																			  *
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
#include "MainWindow.h"
#include "Game.h"
#include <cassert>

Game::Game( MainWindow& wnd )
	:
	wnd( wnd ),
	gfx( wnd )
{
	//Using an initial scalar field to test
	for (int j = 0; j < N; j++)
	{
		for (int i = 0; i < N; i++)
		{
			const float x = float(i - int(N / 2));
			const float y = float(j - int(N / 2));

			//e^(-r^2)
			density[GetId(i, j)] = pow(2.7186f, -0.1f * Vec2(x, y).GetLengthSq());
		}
	}

	//Using an initial vector field to test
	for (int j = 0; j < N; j++)
	{
		for (int i = 0; i < N; i++)
		{
			const float x = float(i - int(N / 2));
			const float y = float(j - int(N / 2));

			//Arbitrary 2D function
			velocity[GetId(i, j)] = Vec2(0.5f*x, y);
		}
	}
}

void Game::Go()
{
	gfx.BeginFrame();	
	UpdateModel();
	ComposeFrame();
	gfx.EndFrame();
}

void Game::UpdateModel()
{
	//Delta time
	const float DT = ft.Mark();
	//Delta time


	//Density source
	if (wnd.mouse.LeftIsPressed())
	{
		Vec2 mousePos(float(wnd.mouse.GetPosX()), float(wnd.mouse.GetPosY()));
		mousePos /= cellDimension;
		density[GetId(int(mousePos.x), int(mousePos.y))] += DT*0.5f;
		if (density[GetId(int(mousePos.x), int(mousePos.y))] > 1.0f)
		{
			density[GetId(int(mousePos.x), int(mousePos.y))] = 1.0f;
		}
	}

	//Add forces
	//Diffuse
	//Diffusion();
	//Move
}

int Game::GetId(int i, int j)
{
	return j * N + i;
}

void Game::DrawDensity()
{
	//Drawing every cell
	for (int j = 0; j < N; j++)
	{
		for (int i = 0; i < N; i++)
		{
			const int A = int(j * cellDimension);
			const int B = int(i * cellDimension);
			for (int J = A; J < A+cellDimension; J++)
			{
				for (int I = B; I < B+cellDimension; I++)
				{
					const int id = GetId(i, j);
					gfx.PutPixel(I, J, Color( char(255 * density[id]), char(255 * density[id]), char(255 * density[id])));
				}
			}
		}
	}
}

void Game::DrawVelocities(bool separated)
{
	//Given that all velocities are normalized we can use this function to draw the vectors
	const float maxDrawnLength = cellDimension / 2.0f;
	if (separated)
	{
		//Finding the max u and v values
		float maxU = velocity[0].x;
		float maxV = velocity[0].y;
		for (int i = 1; i < N * N; i++)
		{
			maxU = std::max<float>(velocity[i].x, maxU);
			maxV = std::max<float>(velocity[i].y, maxV);
		}

		for (int j = 0; j < N; j++)
		{
			for (int i = 0; i < N; i++)
			{
				const Vec2 c(cellDimension * i + maxDrawnLength, cellDimension * j + maxDrawnLength);
				//Horizontal lines
				gfx.DrawLine(c, c + Vec2(velocity[GetId(i,j)].x, 0.0f) * maxDrawnLength / maxU, Colors::Green);
				//Vertical lines
				gfx.DrawLine(c, c + Vec2(0.0f, velocity[GetId(i, j)].y) * maxDrawnLength / maxV, Colors::Red);
			}
		}
	}
	else
	{
		for (int j = 0; j < N; j++)
		{
			for (int i = 0; i < N; i++)
			{
				const Vec2 c(cellDimension * i + maxDrawnLength, cellDimension * j + maxDrawnLength);
				gfx.DrawLine(c, c + velocity[GetId(i, j)].GetNormalized() * maxDrawnLength, Colors::White);
			}
		}
	}
}

void Game::Diffusion()
{
	for (int j = 1; j <= n; j++)
	{
		for (int i = 1; i <= n; i++)
		{
			const int id = GetId(i, j);
			if (densityFlag)
			{
				density[id] = prev_density[id] + diffusionRate *
					(prev_density[id - 1] + prev_density[id + 1] + prev_density[id + N] + prev_density[id - N] - 4 * prev_density[id]);
			}
			else
			{
				prev_density[id] = density[id] + diffusionRate *
					(density[id - 1] + density[id + 1] + density[id + N] + density[id - N] - 4 * density[id]);
			}
		}
	}
}

void Game::ComposeFrame()
{
	DrawDensity();
	//DrawVelocities(true);
}
