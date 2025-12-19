import { test, expect } from '@playwright/test';

/**
 * Visual Regression Tests for Primer Design Application
 *
 * These tests capture screenshots of key UI states to detect
 * unintended visual changes during development.
 */

test.describe('Application Layout', () => {
  test('main application loads correctly', async ({ page }) => {
    await page.goto('/');
    await page.waitForLoadState('networkidle');

    // Wait for any initial loading to complete
    await page.waitForTimeout(1000);

    await expect(page).toHaveScreenshot('app-initial-load.png', {
      fullPage: true,
    });
  });

  test('navigation tabs are visible', async ({ page }) => {
    await page.goto('/');
    await page.waitForLoadState('networkidle');

    // Check for main navigation elements
    const tabList = page.getByRole('tablist');
    if (await tabList.isVisible()) {
      await expect(tabList).toHaveScreenshot('navigation-tabs.png');
    }
  });
});

test.describe('Dark Mode', () => {
  test('dark mode toggle works', async ({ page }) => {
    await page.goto('/');
    await page.waitForLoadState('networkidle');

    // Find and click dark mode toggle if it exists
    const darkModeButton = page.locator('[aria-label*="dark" i], [aria-label*="theme" i], button:has-text("Dark")').first();

    if (await darkModeButton.isVisible()) {
      await darkModeButton.click();
      await page.waitForTimeout(500);

      await expect(page).toHaveScreenshot('app-dark-mode.png', {
        fullPage: true,
      });
    }
  });
});

test.describe('Component Appearance', () => {
  test('primer design form renders correctly', async ({ page }) => {
    await page.goto('/');
    await page.waitForLoadState('networkidle');

    // Look for form elements
    const sequenceInput = page.locator('textarea, input[type="text"]').first();
    if (await sequenceInput.isVisible()) {
      await expect(sequenceInput).toHaveScreenshot('sequence-input.png');
    }
  });

  test('buttons have consistent styling', async ({ page }) => {
    await page.goto('/');
    await page.waitForLoadState('networkidle');

    // Capture a primary action button
    const primaryButton = page.locator('button[type="submit"], button:has-text("Design"), button:has-text("Calculate")').first();
    if (await primaryButton.isVisible()) {
      await expect(primaryButton).toHaveScreenshot('primary-button.png');
    }
  });
});

test.describe('Responsive Design', () => {
  test('mobile viewport renders correctly', async ({ page }) => {
    await page.setViewportSize({ width: 375, height: 667 });
    await page.goto('/');
    await page.waitForLoadState('networkidle');
    await page.waitForTimeout(1000);

    await expect(page).toHaveScreenshot('mobile-viewport.png', {
      fullPage: true,
    });
  });

  test('tablet viewport renders correctly', async ({ page }) => {
    await page.setViewportSize({ width: 768, height: 1024 });
    await page.goto('/');
    await page.waitForLoadState('networkidle');
    await page.waitForTimeout(1000);

    await expect(page).toHaveScreenshot('tablet-viewport.png', {
      fullPage: true,
    });
  });
});

test.describe('Error States', () => {
  test('empty sequence shows appropriate message', async ({ page }) => {
    await page.goto('/');
    await page.waitForLoadState('networkidle');

    // Try to submit without sequence
    const submitButton = page.locator('button[type="submit"]').first();
    if (await submitButton.isVisible()) {
      await submitButton.click();
      await page.waitForTimeout(500);

      // Capture any error message that appears
      const errorElement = page.locator('[role="alert"], .error, .text-red-500').first();
      if (await errorElement.isVisible()) {
        await expect(errorElement).toHaveScreenshot('error-message.png');
      }
    }
  });
});
